use rusqlite::{params, Connection, Result};
use std::io::{BufRead, BufReader, Read, Write};
use std::sync::{mpsc, Arc, Mutex};
use std::thread;
use std::io;
use std::path::{Path, PathBuf};
use log::{debug, info};
use std::fs::{self, File};
use flate2::bufread::GzDecoder;
use rayon::{prelude::*, ThreadPoolBuilder};
use std::time::Instant;

fn open_fastq_file(fastq_file: &Path) -> Result<BufReader<Box<dyn Read>>, std::io::Error> {
    if fastq_file.extension().map(|ext| ext == "gz").unwrap_or(false) {
        let file = File::open(fastq_file)?;
        let decoder = GzDecoder::new(BufReader::new(file));
        Ok(BufReader::new(Box::new(decoder) as Box<dyn Read>))
    } else {
        let file = File::open(fastq_file)?;
        Ok(BufReader::new(Box::new(file) as Box<dyn Read>))
    }
}

fn to_io_error(e: rusqlite::Error) -> io::Error {
    io::Error::new(io::ErrorKind::Other, e.to_string())
}

pub fn indexfastqreads(readdir: &PathBuf, outputpath: &PathBuf, threads: usize) -> io::Result<()> {
    let start_time = Instant::now();
    let db_path = outputpath.join("reads.db");

    debug!("db_path {:?}", db_path);
    
    if db_path.exists() {
        info!("Database already exists at: {}", db_path.display());
    } else {
        info!("Creating a read index database: {}", db_path.display());
    }
    
    let mut fastq_files = Vec::new();
    if let Ok(entries) = fs::read_dir(&readdir) {
        for entry in entries.filter_map(Result::ok) {
            let path = entry.path();
            if path.extension()
                .map(|ext| ext == "fastq" || ext == "gz")
                .unwrap_or(false)
            {
                fastq_files.push(path);
            }
        }
    }

    if let Some(first) = fastq_files.get(0) {
        debug!("The first element is: {:?}", first);
    } else {
        debug!("The vector is empty.");
    }

    // Create the table if it doesn't exist
    if !db_path.exists() {
        let conn = Arc::new(Mutex::new(Connection::open(&db_path).map_err(to_io_error)?));
        {
            debug!("insider conn setting");
            let conn = conn.try_lock().map_err(|e| {
                let err_msg = format!("Failed to lock connection: {}", e);
                io::Error::new(io::ErrorKind::Other, err_msg)
            })?;
            std::io::stdout().flush().unwrap();
            debug!("initialized conn try lock");
            match conn.prepare("PRAGMA journal_mode = WAL;") {
                Ok(mut stmt) => {
                    let result = stmt.query([]); // PRAGMA statements that return results
                    match result {
                        Ok(mut rows) => {
                            if let Ok(Some(row)) = rows.next() {
                                // Consume the row (you can extract the result if needed)
                                debug!("PRAGMA journal_mode executed successfully: {:?}", row);
                            }
                        }
                        Err(e) => {
                            debug!("Error executing PRAGMA journal_mode: {}", e);
                            return Err(io::Error::new(io::ErrorKind::Other, e.to_string()));
                        }
                    }
                }
                Err(e) => {
                    debug!("Error preparing PRAGMA journal_mode: {}", e);
                    return Err(io::Error::new(io::ErrorKind::Other, e.to_string()));
                }
            }
            debug!("pragma journal mode");
            
            if let Err(e) = conn.execute("PRAGMA synchronous = NORMAL;", []) {
                debug!("Error executing PRAGMA synchronous: {}", e);
                return Err(io::Error::new(io::ErrorKind::Other, e.to_string()));
            }
            debug!("pragma sync");
            conn.execute("PRAGMA temp_store = MEMORY;", []).map_err(to_io_error)?;
            debug!("pragma temp store");
            
            match conn.prepare("PRAGMA mmap_size = 10485760;") {
                Ok(mut stmt) => {
                    let result = stmt.query([]); // PRAGMA statements that return results
                    match result {
                        Ok(mut rows) => {
                            // Consume the row, if needed
                            if let Ok(Some(row)) = rows.next() {
                                debug!("PRAGMA mmap_size executed successfully: {:?}", row);
                            }
                        }
                        Err(e) => {
                            eprintln!("Error executing PRAGMA mmap_size: {}", e);
                            return Err(io::Error::new(io::ErrorKind::Other, e.to_string()));
                        }
                    }
                }
                Err(e) => {
                    eprintln!("Error preparing PRAGMA mmap_size: {}", e);
                    return Err(io::Error::new(io::ErrorKind::Other, e.to_string()));
                }
            }

            debug!("pragma mmap size 10485760");
            conn.execute(
                "CREATE TABLE IF NOT EXISTS fastq (
                    file_id TEXT NOT NULL,
                    read_id TEXT NOT NULL,
                    sequence TEXT NOT NULL,
                    quality TEXT NOT NULL,
                    PRIMARY KEY (file_id, read_id)
                    );",
                    [],
                ).map_err(to_io_error)?;
            debug!("conn execute");
            }
        debug!("after conn setting");    
        // Create the channel for communication between threads
        let (tx, rx) = mpsc::sync_channel::<(String, String, String, String)>(100);
        debug!("tx and rx is initialized");
        let conn_writer = Arc::clone(&conn);
        debug!("conn writer initialized");
        let writer_handle = thread::spawn(move || {
            let mut conn = conn_writer.lock().unwrap();
            debug!("mut conn inside writer handle");
            let transaction = conn.transaction().unwrap();
            {

                debug!("transaction unwrap");
                let mut insert_stmt = transaction.prepare(
                    "INSERT OR IGNORE INTO fastq (file_id, read_id, sequence, quality) VALUES (?, ?, ?, ?)",
                ).unwrap();
                debug!("insert stmt setup");
                let mut batch: Vec<(String, String, String, String)> = Vec::new();
                debug!("batch is created");
                while let Ok(item) = rx.recv() {
                    batch.push(item);
                    if batch.len() >= 10000 {
                        for (file_id, read_id, sequence, quality) in batch.drain(..) {
                            if let Err(e) = insert_stmt.execute(params![file_id, read_id, sequence, quality]) {
                                eprintln!("Error inserting data: {:?}", e);
                            }
                        }
                    }
                }
                debug!("batch is complete after while loop");
                // Process any remaining items
                for (file_id, read_id, sequence, quality) in batch {
                    if let Err(e) = insert_stmt.execute(params![file_id, read_id, sequence, quality]) {
                        eprintln!("Error inserting data: {:?}", e);
                    }
                }
            }
            transaction.commit().unwrap();
            debug!("Transaction is completed");
        });
      
        ThreadPoolBuilder::new()
            .num_threads(fastq_files.len().min(threads))
            .build_global()
            .expect("Failed to create a thread pool");

        fastq_files.par_iter().for_each(|fastq_file| {
            if let Err(err) = process_fastq_file(fastq_file, &tx) {
                eprintln!("Error processing file {:?}: {}", fastq_file, err);
            }
        });
        
        drop(tx); // Close the sender channel.
        writer_handle.join().expect("Writer thread panicked");
        
        let conn = conn.lock().map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
        let count: usize = conn.query_row("SELECT COUNT(*) FROM fastq", [], |row| row.get(0))
            .map_err(to_io_error)?;
        debug!("Rows in fastq table: {}", count);
    }
    let duration = start_time.elapsed();
    info!("Time taken for the function: {:?}", duration);
    Ok(())
}

fn process_fastq_file(
    fastq_file: &PathBuf,
    tx: &mpsc::SyncSender<(String, String, String, String)>,
) -> io::Result<()> {
    let fastq_reader = open_fastq_file(fastq_file)?;
    let mut lines = fastq_reader.lines();
    let mut batch = Vec::with_capacity(1000);
    debug!("Processing fastq file");
    while let (Some(Ok(header1)), Some(Ok(seq1)), Some(Ok(_)), Some(Ok(qual1))) =
        (lines.next(), lines.next(), lines.next(), lines.next())
    {
        batch.push((
            fastq_file.to_string_lossy().to_string(),
            header1.trim_start_matches('@').to_string(),
            seq1,
            qual1,
        ));

        for record in batch.drain(..) {
            let _ = tx.send(record);
        }
    }

    Ok(())
}


// // Spawn worker threads to process FASTQ files
// fastq_files.par_iter().for_each(|fastq_file| {
//     // let file = File::open(fastq_file).unwrap();
//     // let mut decoder = GzDecoder::new(BufReader::new(file));
//     // let mut decompress_file = Vec::new();
//     // decoder.read_to_end(&mut decompress_file).unwrap();
//     // let fastq_reader = BufReader::new(decompress_file.as_slice());
//     let result: Result<(), std::io::Error> = (|| {
//         let fastq_reader = open_fastq_file(fastq_file)?;

//         let mut lines = fastq_reader.lines();
//         while let (Some(Ok(header1)), Some(Ok(seq1)), Some(Ok(_)), Some(Ok(qual1))) = (
//             lines.next(),
//             lines.next(),
//             lines.next(),
//             lines.next(),
//         ) {
//             tx.send((
//                 fastq_file.to_string_lossy().to_string(),
//                 header1.trim_start_matches('@').to_string(),
//                 seq1,
//                 qual1,
//             ))
//             .unwrap();
//         }
//         Ok(())
//     })();
    
//     // Handle the result of the closure (e.g., log error or handle it)
//     if let Err(err) = result {
//         eprintln!("Error processing file {:?}: {}", fastq_file, err);
//     }
// });