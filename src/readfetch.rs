use std::arch::x86_64::_CMP_NEQ_US;
use std::collections::HashSet;
use std::process::{Command, Stdio};
// use clap::error;
// use flate2::bufread::GzDecoder;
use hashbrown::HashMap;
use std::path::PathBuf;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, Write};
use rusqlite::Result;
// use std::io::{self, BufRead, BufReader, Write, Read};
// use rusqlite::{params, Connection, Result};
use crate::utility;
use log::debug;

fn read_mapid_file(mapid_file: &str) -> io::Result<HashMap<String, Vec<String>>> {
    let mut scaffold_to_reads: HashMap<String, Vec<String>> = HashMap::new();
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
            scaffold_to_reads
                .entry(scaffold_id.to_owned())
                .or_insert_with(Vec::new)
                .push(read_id.to_owned());
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

    Ok(scaffold_to_reads)

}


fn write_selected_reads(
    fastq_files: Vec<String>,
    enriched_scaffolds: &HashSet<String>,
    mapid_file: &str,
    output_fastq: &Vec<String>,
    is_paired: bool,
    create_new: bool,
) -> Result<(), io::Error> {
    // Extract reads mapped to enriched scaffolds
    let scaffold_to_reads = read_mapid_file(mapid_file)?;
    
    let selected_contigs: HashSet<String> = scaffold_to_reads
        .clone()
        .into_iter()
        .filter(|(scaffold, _)| enriched_scaffolds.contains(scaffold))
        .map(|(contig,_)| contig)
        .collect();
    debug!("selected contigs: {}", selected_contigs.len());

    let selected_reads: HashSet<String> = scaffold_to_reads
        .into_iter()
        .filter(|(scaffold, _)| enriched_scaffolds.contains(scaffold))
        .flat_map(|(_, read)| read)
        .collect();

    // Prepare the read ID file
    let readid_file = format!("{}_readid", output_fastq[0].replace(".fastq", ""));
    debug!(
        "Number of selected reads: {:?}, Read ID file: {:?}",
        selected_reads.len(),
        readid_file
    );
    
    let mut idfile = if create_new {
        File::create(&readid_file)?
    } else {
        OpenOptions::new().create(true).append(true).open(&readid_file)?
    };

    for read in &selected_reads {
        writeln!(idfile, "{}", read)?;
    }

    // Define the file writing logic based on `create_new`
    let write_to_file = |file_path: &String, data: &[u8]| -> Result<(), io::Error> {
        let mut file = if create_new {
            File::create(file_path)?
        } else {
            OpenOptions::new().create(true).append(true).open(file_path)?
        };
        file.write_all(data)
    };

    if is_paired {
        // Process forward reads
        let output1 = Command::new("seqtk")
            .arg("subseq")
            .arg(&fastq_files[0])
            .arg(&readid_file)
            .stdout(Stdio::piped())
            .output()
            .expect("Failed to execute seqtk for forward reads");
        write_to_file(&output_fastq[0], &output1.stdout)?;

        // Process reverse reads
        let output2 = Command::new("seqtk")
            .arg("subseq")
            .arg(&fastq_files[1])
            .arg(&readid_file)
            .stdout(Stdio::piped())
            .output()
            .expect("Failed to execute seqtk for reverse reads");
        write_to_file(&output_fastq[1], &output2.stdout)?;
    } else {
        // Process single-end reads
        let output = Command::new("seqtk")
            .arg("subseq")
            .arg(&fastq_files[0])
            .arg(&readid_file)
            .stdout(Stdio::piped())
            .output()
            .expect("Failed to execute seqtk for single-end reads");
        write_to_file(&output_fastq[0], &output.stdout)?;
    }

    // if let Some(first_element) = selected_reads.iter().next() {
    //     debug!("The first element is: {}", first_element);
    // } else {
    //     debug!("The HashSet is empty.");
    // }

    // // Prepare the output file
    // let output_file = if create_new {
    //     File::create(&output_fastq)? // Create a new file if flag is true
    // } else {
    //     OpenOptions::new()
    //         .create(true)
    //         .append(true)  // Open in append mode if flag is false
    //         .open(&output_fastq)?
    // };

    // let mut output_file = io::BufWriter::new(output_file);
    // let conn = Connection::open(db_path).map_err(|e| {
    //     io::Error::new(io::ErrorKind::Other, format!("Database connection failed: {}", e))
    // })?;

    // for read_id in selected_reads {
    //     let mut query = String::from("SELECT file_id, read_id, sequence, quality FROM fastq WHERE read_id LIKE ?1 || '.%'");
    //     let mut params: Vec<&dyn rusqlite::ToSql> = vec![&read_id];
    
    //     if !fastq_files.is_empty() {
    //         let placeholders: String = fastq_files
    //             .iter()
    //             .enumerate()
    //             .map(|(i, _)| format!("?{}", i + 2))
    //             .collect::<Vec<_>>()
    //             .join(", ");
            
    //         query.push_str(&format!(" AND file_id IN ({})", placeholders));
        
    //         // Add file names to parameters
    //         for file_name in &fastq_files {
    //             params.push(file_name);
    //         }
    //     }
    
    //     // Prepare statement and bind parameters
    //     let mut stmt = conn.prepare(&query).map_err(|e| {
    //         io::Error::new(io::ErrorKind::Other, format!("Query preparation failed: {}", e))
    //     })?;
    
    //     // Execute the query with the prepared parameters
    //     let mut rows = stmt.query(params.as_slice()).map_err(|e| {
    //         io::Error::new(io::ErrorKind::Other, format!("Query execution failed: {}", e))
    //     })?;
    
    //     // Write results to the file
    //     while let Some(row) = rows.next().map_err(|e| {
    //         io::Error::new(io::ErrorKind::Other, format!("Row retrieval failed: {}", e))
    //     })? {
    //         let read_id: String = row.get(1)
    //             .map_err(|e|
    //             io::Error::new(
    //             io::ErrorKind::Other,
    //             e.to_string()))?;
    //         let sequence: String = row.get(2)
    //             .map_err(|e|
    //             io::Error::new(
    //             io::ErrorKind::Other,
    //             e.to_string()))?;
    //         let quality: String = row.get(3)
    //             .map_err(|e|
    //             io::Error::new(
    //             io::ErrorKind::Other,
    //             e.to_string()))?;
    
    //         writeln!(
    //             output_file,
    //             "@{}\n{}\n+\n{}", read_id, sequence, quality
    //         )?;
    //     }
    // }
    Ok(())
}

pub fn fetch_fastqreads(
    enriched_scaffolds: &HashSet<String>,
    mapids: &str,
    fastq_files: Vec<String>,
    outputbin: PathBuf,
    is_paired: bool,
    create_new: bool,
) -> Result<(), Box<dyn std::error::Error>> {

    let output_fastq: Vec<String> = if is_paired {
            let read_file1 = PathBuf::from(format!("{}",
                utility::get_output_binname(outputbin.to_str().expect(""))
                .to_str()
                .expect("Invalid UTF-8 in file path")
                .replace(".fasta", "_1.fastq")));
            let read_file2 = PathBuf::from(format!("{}",
                utility::get_output_binname(outputbin.to_str().expect(""))
                .to_str()
                .expect("Invalid UTF-8 in file path")
                .replace(".fasta", "_2.fastq")));
    
            vec![
                read_file1.to_str().expect("Failed to convert PathBuf to &str").to_string(),
                read_file2.to_str().expect("Failed to convert PathBuf to &str").to_string()
            ]
        } else {
            let read_file = PathBuf::from(format!("{}",
                utility::get_output_binname(outputbin.to_str().expect(""))
                .to_str()
                .expect("Invalid UTF-8 in file path")
                .replace(".fasta", ".fastq")));
            
            vec![
                read_file.to_str().expect("Failed to convert PathBuf to &str").to_string()
            ]
        };

    debug!("Writing reads for {:?}", output_fastq);
    let _ = write_selected_reads(
        fastq_files,
        &enriched_scaffolds,
        mapids,
        &output_fastq,
        is_paired,
        create_new);

    let file_metadata = std::fs::metadata(&output_fastq[0]).map_err(|e| {
        eprintln!("Error accessing {:?}: {}", output_fastq, e);
        Box::new(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("Failed to access {:?}", output_fastq),
        )) as Box<dyn std::error::Error>
    })?;

    if file_metadata.len() == 0 {
        eprintln!("Error: The output file {:?} is empty", output_fastq[0]);
        return Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("No reads were written to {:?}", output_fastq[0]),
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