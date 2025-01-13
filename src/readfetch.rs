use std::collections::HashSet;
// use clap::error;
// use flate2::bufread::GzDecoder;
use hashbrown::HashMap;
use std::path::PathBuf;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, Write};
use rusqlite::{Connection, Result};
// use std::io::{self, BufRead, BufReader, Write, Read};
// use rusqlite::{params, Connection, Result};
use crate::utility;
use log::debug;

fn read_mapid_file(mapid_file: &str) -> io::Result<HashMap<String, String>> {
    let mut read_to_scaffold = HashMap::new();
    let file = File::open(mapid_file).map_err(|e| {
        io::Error::new(io::ErrorKind::NotFound,
            format!("Failed to open file '{}': {}", mapid_file, e))
    })?;

    let reader = io::BufReader::new(file);

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| {
            io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Failed to read line {} in '{}': {}", line_num + 1, mapid_file, e),
            )
        })?;

        let mut parts = line.split_whitespace();
        let read_id = parts.next();
        let scaffold_id = parts.next();

        if let (Some(read_id), Some(scaffold_id)) = (read_id, scaffold_id) {
            // Insert directly into hashbrown's HashMap
            read_to_scaffold.insert(read_id.to_owned(), scaffold_id.to_owned());
        } else {
            return Err(
                io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Invalid format at line {} in '{}': '{}'",
                line_num + 1,
                mapid_file,
                line),
            ));
        }
    }

    Ok(read_to_scaffold)

}

// fn write_selected_reads(
//     fastq_file: &str,
//     fastq_file2: Option<&str>,
//     enriched_scaffolds: &HashSet<String>,
//     mapid_file: &str,
//     output_fastq: PathBuf,
//     is_paired: bool,
//     create_new: bool,
// ) -> io::Result<()> {
//     let read_to_scaffold = read_mapid_file(mapid_file)?;
//     let selected_reads: HashSet<String> = read_to_scaffold
//         .into_iter()
//         .filter(|(_, scaffold)| enriched_scaffolds.contains(scaffold))
//         .map(|(read, _)| read)
//         .collect();
    
//     debug!("length of selected reads {:?}", selected_reads.len());

//     if let Some(first_element) = selected_reads.iter().next() {
//         debug!("The first element is: {}", first_element);
//     } else {
//         debug!("The HashSet is empty.");
//     }

//     // Open the first fastq file
//     let file = File::open(fastq_file).map_err(|e| {
//         io::Error::new(io::ErrorKind::NotFound,
//             format!("Failed to open file '{}': {}", fastq_file, e))
//     })?;

//     let mut decoder = GzDecoder::new(BufReader::new(file));
//     let mut decompress_file = Vec::with_capacity(1_000_000_000);

//     decoder.read_to_end(&mut decompress_file)?;

//     let fastq_reader = BufReader::new(decompress_file.as_slice());

//     // Prepare the output file
//     let output_file = if create_new {
//         File::create(&output_fastq)? // Create a new file if flag is true
//     } else {
//         OpenOptions::new()
//             .create(true)
//             .append(true)  // Open in append mode if flag is false
//             .open(&output_fastq)?
//     };

//     let mut output_file = io::BufWriter::new(output_file);

//     if is_paired {
//         // Handle paired-end reads from a second file
//         let file2 = File::open(fastq_file2.expect("Paired-end flag set but no second file provided"))?;
       
//         let mut decoder2 = GzDecoder::new(BufReader::new(file2));
//         let mut decompress_file2 = Vec::with_capacity(1_000_000_000);

//         decoder2.read_to_end(&mut decompress_file2)?;

//         let fastq_reader2 = BufReader::new(decompress_file2.as_slice());

//         let mut reader1 = fastq_reader.lines();
//         let mut reader2 = fastq_reader2.lines();

//         loop {
//             let header1 = match reader1.next() {
//                 Some(Ok(line)) => line,
//                 Some(Err(e)) => return Err(
//                     io::Error::new(
//                     io::ErrorKind::Other, format!("Error reading FASTQ file 1: {}", e))),
//                 None => break, // EOF
//             };
    
//             let seq1 = reader1
//                 .next()
//                 .unwrap_or_else(|| 
//                     Err(
//                     io::Error::new(
//                     io::ErrorKind::UnexpectedEof, "Unexpected EOF in FASTQ file 1")))?;
//             let plus1 = reader1
//                 .next()
//                 .unwrap_or_else(||
//                     Err(
//                     io::Error::new(
//                     io::ErrorKind::UnexpectedEof, "Unexpected EOF in FASTQ file 1")))?;
//             let qual1 = reader1
//                 .next()
//                 .unwrap_or_else(||
//                     Err(
//                     io::Error::new(
//                     io::ErrorKind::UnexpectedEof, "Unexpected EOF in FASTQ file 1")))?;
    
//             let header2 = match reader2.next() {
//                 Some(Ok(line)) => line,
//                 Some(Err(e)) => return Err(
//                     io::Error::new(
//                     io::ErrorKind::Other, format!("Error reading FASTQ file 2: {}", e))),
//                 None => break, // EOF
//             };
    
//             let seq2 = reader2
//                 .next()
//                 .unwrap_or_else(||
//                 Err(
//                 io::Error::new(
//                 io::ErrorKind::UnexpectedEof, "Unexpected EOF in FASTQ file 2")))?;
//             let plus2 = reader2
//                 .next()
//                 .unwrap_or_else(||
//                 Err(
//                 io::Error::new(
//                 io::ErrorKind::UnexpectedEof, "Unexpected EOF in FASTQ file 2")))?;
//             let qual2 = reader2
//                 .next()
//                 .unwrap_or_else(||
//                 Err(
//                 io::Error::new(
//                 io::ErrorKind::UnexpectedEof, "Unexpected EOF in FASTQ file 2")))?;
    
//             // Process the trimmed read ID
//             let read_id_trimmed = header1
//                 .trim_start_matches('@')
//                 .trim_end_matches("/1")
//                 .to_string();
    
//             // Check if the read is in the selected set and write if matched
//             if selected_reads.contains(&read_id_trimmed) {
//                 writeln!(output_file, "{}/1", header1)?;
//                 writeln!(output_file, "{}", seq1)?;
//                 writeln!(output_file, "{}", plus1)?;
//                 writeln!(output_file, "{}", qual1)?;
    
//                 writeln!(output_file, "{}/2", header2)?;
//                 writeln!(output_file, "{}", seq2)?;
//                 writeln!(output_file, "{}", plus2)?;
//                 writeln!(output_file, "{}", qual2)?;
//             }
//         }
//     } else {
//         let mut reader = fastq_reader.lines();
    
//         loop {
//             // Read 4 lines per record
//             let header = match reader.next() {
//                 Some(Ok(line)) => line,
//                 Some(Err(e)) => return Err(
//                     io::Error::new(
//                     io::ErrorKind::Other, format!("Error reading FASTQ file: {}", e))),
//                 None => break, // EOF
//             };
    
//             let seq = reader
//                 .next()
//                 .unwrap_or_else(||
//                 Err(
//                 io::Error::new(
//                 io::ErrorKind::UnexpectedEof, "Unexpected EOF in FASTQ file")))?;
//             let plus = reader
//                 .next()
//                 .unwrap_or_else(||
//                 Err(
//                 io::Error::new(
//                 io::ErrorKind::UnexpectedEof, "Unexpected EOF in FASTQ file")))?;
//             let qual = reader
//                 .next()
//                 .unwrap_or_else(||
//                 Err(
//                 io::Error::new(
//                 io::ErrorKind::UnexpectedEof, "Unexpected EOF in FASTQ file")))?;
    
//             // Process the trimmed read ID
//             let read_id_trimmed = header
//                 .trim_start_matches('@')
//                 .trim_end_matches("/1")
//                 .trim_end_matches("/2")
//                 .to_string();
    
//             // Check if the read is in the selected set and write if matched
//             if selected_reads.contains(&read_id_trimmed) {
//                 writeln!(output_file, "{}", header)?;
//                 writeln!(output_file, "{}", seq)?;
//                 writeln!(output_file, "{}", plus)?;
//                 writeln!(output_file, "{}", qual)?;
//             }
//         }
//     }
//     output_file.flush()?;
//     Ok(())
// }


fn write_selected_reads(
    fastq_files: Vec<String>,
    db_path: &PathBuf,
    enriched_scaffolds: &HashSet<String>,
    mapid_file: &str,
    output_fastq: PathBuf,
    create_new: bool,
) -> Result<(), std::io::Error> {
    let read_to_scaffold = read_mapid_file(mapid_file)?;
    let selected_reads: HashSet<String> = read_to_scaffold
        .into_iter()
        .filter(|(_, scaffold)| enriched_scaffolds.contains(scaffold))
        .map(|(read, _)| read)
        .collect();
    
    debug!("length of selected reads {:?}", selected_reads.len());

    if let Some(first_element) = selected_reads.iter().next() {
        debug!("The first element is: {}", first_element);
    } else {
        debug!("The HashSet is empty.");
    }

    // Prepare the output file
    let output_file = if create_new {
        File::create(&output_fastq)? // Create a new file if flag is true
    } else {
        OpenOptions::new()
            .create(true)
            .append(true)  // Open in append mode if flag is false
            .open(&output_fastq)?
    };

    let mut output_file = io::BufWriter::new(output_file);
    let conn = Connection::open(db_path).map_err(|e| {
        io::Error::new(io::ErrorKind::Other, format!("Database connection failed: {}", e))
    })?;

    for read_id in selected_reads {
        let mut query = String::from("SELECT file_id, read_id, sequence, quality FROM fastq WHERE read_id LIKE ?1 || '.%'");
        let mut params: Vec<&dyn rusqlite::ToSql> = vec![&read_id];
    
        if !fastq_files.is_empty() {
            let placeholders: String = fastq_files
                .iter()
                .enumerate()
                .map(|(i, _)| format!("?{}", i + 2))
                .collect::<Vec<_>>()
                .join(", ");
            
            query.push_str(&format!(" AND file_id IN ({})", placeholders));
        
            // Add file names to parameters
            for file_name in &fastq_files {
                params.push(file_name);
            }
        }
    
        // Prepare statement and bind parameters
        let mut stmt = conn.prepare(&query).map_err(|e| {
            io::Error::new(io::ErrorKind::Other, format!("Query preparation failed: {}", e))
        })?;
    
        // Execute the query with the prepared parameters
        let mut rows = stmt.query(params.as_slice()).map_err(|e| {
            io::Error::new(io::ErrorKind::Other, format!("Query execution failed: {}", e))
        })?;
    
        // Write results to the file
        while let Some(row) = rows.next().map_err(|e| {
            io::Error::new(io::ErrorKind::Other, format!("Row retrieval failed: {}", e))
        })? {
            let read_id: String = row.get(1)
                .map_err(|e|
                io::Error::new(
                io::ErrorKind::Other,
                e.to_string()))?;
            let sequence: String = row.get(2)
                .map_err(|e|
                io::Error::new(
                io::ErrorKind::Other,
                e.to_string()))?;
            let quality: String = row.get(3)
                .map_err(|e|
                io::Error::new(
                io::ErrorKind::Other,
                e.to_string()))?;
    
            writeln!(
                output_file,
                "@{}\n{}\n+\n{}", read_id, sequence, quality
            )?;
        }
    }
    Ok(())
}

pub fn fetch_fastqreads(
    enriched_scaffolds: &HashSet<String>,
    mapids: &str,
    fastq_files: Vec<String>,
    db_path: &PathBuf,
    outputbin: PathBuf,
    create_new: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let output_fastq = PathBuf::from(format!("{}",
        utility::get_output_binname(outputbin.to_str().expect(""))
        .to_str()
        .expect("Invalid UTF-8 in file path")
        .replace(".fasta", ".fastq")));
    debug!("Writing reads for {:?}", output_fastq);
    let _ = write_selected_reads(
        fastq_files,
        db_path,
        &enriched_scaffolds,
        mapids,
        output_fastq.clone(),
        create_new);

    let file_metadata = std::fs::metadata(&output_fastq).map_err(|e| {
        eprintln!("Error accessing {:?}: {}", output_fastq, e);
        Box::new(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("Failed to access {:?}", output_fastq),
        )) as Box<dyn std::error::Error>
    })?;

    if file_metadata.len() == 0 {
        eprintln!("Error: The output file {:?} is empty", output_fastq);
        return Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("No reads were written to {:?}", output_fastq),
        )));
    }
    Ok(())
}

// pub fn fetch_fastqreads(
//     enriched_scaffolds: &HashSet<String>,
//     mapids: &str,
//     read_fastq: &str,
//     read_fastq2: Option<&str>,
//     outputbin: PathBuf,
//     is_paired: bool,
//     create_new: bool,
// ) {
//     let output_fastq = PathBuf::from(format!("{}",
//         utility::get_output_binname(outputbin.to_str().expect(""))
//         .to_str()
//         .expect("Invalid UTF-8 in file path")
//         .replace(".fasta", ".fastq")));
//     debug!("Writing reads for {:?}", output_fastq);
//     let _ = write_selected_reads(
//         &read_fastq,
//         read_fastq2,
//         &enriched_scaffolds,
//         mapids,
//         output_fastq,
//         is_paired,
//         create_new);
// }