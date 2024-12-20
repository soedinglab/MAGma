use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, Write};
use flate2::read::GzDecoder;
use crate::utility;
use log::debug;

fn read_mapid_file(mapid_file: &str) -> io::Result<HashMap<String, String>> {
    let mut read_to_scaffold = HashMap::new();
    let file = File::open(mapid_file).map_err(|e| {
        io::Error::new(io::ErrorKind::NotFound,
            format!("Failed to open file '{}': {}", mapid_file, e))
    })?;

    let reader = io::BufReader::new(file);

    for line in reader.lines() {
        let line = line.map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidData,
            format!("Failed to read line in '{}': {}", mapid_file, e))
        })?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 2 {
            let read_id = parts[0].split('.')
                .next().unwrap_or("").to_string();
            let scaffold_id = parts[1].to_string();
            read_to_scaffold.insert(read_id, scaffold_id);
        } else {
            // Return error if the line doesn't contain exactly two parts
            return Err(io::Error::new(
                io::ErrorKind::InvalidData, format!(
                "Invalid line format in '{}': {}", mapid_file, line
            )));
        }
    }
    Ok(read_to_scaffold)
}

fn write_selected_reads(
    fastq_file: &str,
    fastq_file2: Option<&str>,
    enriched_scaffolds: &HashSet<String>,
    mapid_file: &str,
    output_fastq: PathBuf,
    is_paired: bool,
    create_new: bool,
) -> io::Result<()> {
    let read_to_scaffold = read_mapid_file(mapid_file)?;
    let selected_reads: HashSet<String> = read_to_scaffold
        .iter()
        .filter(|(_, scaffold)| enriched_scaffolds.contains(*scaffold))
        .map(|(read, _)| read.clone())
        .collect();
    let file = File::open(fastq_file).map_err(|e| {
        io::Error::new(io::ErrorKind::NotFound,
            format!("Failed to open file '{}': {}", fastq_file, e))
    })?;
    let reader: Box<dyn BufRead> = if fastq_file.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let output_file = if create_new {
        File::create(&output_fastq)? // Create a new file if flag is true
    } else {
        OpenOptions::new()
            .create(true)
            .append(true)  // Open in append mode if flag is false
            .open(&output_fastq)?
    };

    let mut output_file = io::BufWriter::new(output_file);

    // let mut lines = reader.lines();
    // while let Some(header) = lines.next() {
    //     let sequence = lines.next().unwrap_or(Ok(String::new()))?;
    //     let plus_line = lines.next().unwrap_or(Ok(String::new()))?;
    //     let quality = lines.next().unwrap_or(Ok(String::new()))?;
    //     let read_id = header?.to_string();
    //     let read_id_trimmed = read_id
    //         .trim_start_matches('@')
    //         .trim_end_matches("/1")
    //         .trim_end_matches("/2")
    //         .to_string();
    //     if selected_reads.contains(&read_id_trimmed) {
    //         writeln!(output_file, "{}", read_id)?;
    //         writeln!(output_file, "{}", sequence)?;
    //         writeln!(output_file, "{}", plus_line)?;
    //         writeln!(output_file, "{}", quality)?;
    //     }
    // }
    if is_paired {
        // Handle paired-end reads
        let file2 = File::open(fastq_file2.expect("Paired-end flag set but no second file provided"))?;
        let reader2: Box<dyn BufRead> = if fastq_file2.unwrap().ends_with(".gz") {
            Box::new(BufReader::new(GzDecoder::new(file2)))
        } else {
            Box::new(BufReader::new(file2))
        };

        let mut lines1 = reader.lines();
        let mut lines2 = reader2.lines();

        while let (Some(header1), Some(header2)) = (lines1.next(), lines2.next()) {
            let header1 = header1?;
            let header2 = header2?;
            if header1.starts_with('@') && header2.starts_with('@') {
                let sequence1 = lines1.next()
                    .unwrap_or(Ok(String::new()))?;
                let plus_line1 = lines1.next()
                    .unwrap_or(Ok(String::new()))?;
                let quality1 = lines1.next()
                    .unwrap_or(Ok(String::new()))?;

                let sequence2 = lines2.next()
                    .unwrap_or(Ok(String::new()))?;
                let plus_line2 = lines2.next()
                    .unwrap_or(Ok(String::new()))?;
                let quality2 = lines2.next()
                    .unwrap_or(Ok(String::new()))?;

                let read_id_trimmed = header1
                    .trim_start_matches('@')
                    .trim_end_matches("/1")
                    .to_string();

                if selected_reads.contains(&read_id_trimmed) {
                    writeln!(output_file, "{}", header1)?;
                    writeln!(output_file, "{}", sequence1)?;
                    writeln!(output_file, "{}", plus_line1)?;
                    writeln!(output_file, "{}", quality1)?;

                    writeln!(output_file, "{}", header2)?;
                    writeln!(output_file, "{}", sequence2)?;
                    writeln!(output_file, "{}", plus_line2)?;
                    writeln!(output_file, "{}", quality2)?;
                }
            }
        }
    } else {
        let mut lines = reader.lines();
        while let Some(header) = lines.next() {
            let header = header?;
            if header.starts_with('@') {
                let sequence = lines.next()
                    .unwrap_or(Ok(String::new()))?;
                let plus_line = lines.next()
                    .unwrap_or(Ok(String::new()))?;
                let quality = lines.next()
                    .unwrap_or(Ok(String::new()))?;

                let read_id_trimmed = header
                    .trim_start_matches('@')
                    .trim_end_matches("/1")
                    .trim_end_matches("/2")
                    .to_string();

                if selected_reads.contains(&read_id_trimmed) {
                    writeln!(output_file, "{}", header)?;
                    writeln!(output_file, "{}", sequence)?;
                    writeln!(output_file, "{}", plus_line)?;
                    writeln!(output_file, "{}", quality)?;
                }
            }
        }
    }

    Ok(())
}

pub fn fetch_fastqreads(
    enriched_scaffolds: &HashSet<String>,
    mapids: &str,
    read_fastq: &str,
    read_fastq2: Option<&str>,
    outputbin: PathBuf,
    is_paired: bool,
    create_new: bool,
) {
    let output_fastq = PathBuf::from(format!("{}",
        utility::get_output_scaffoldname(outputbin.to_str().expect(""))
        .to_str()
        .expect("Invalid UTF-8 in file path")
        .replace(".fasta", ".fastq")));
    debug!("Writing reads for {:?}", output_fastq);
    let _ = write_selected_reads(
        &read_fastq,
        read_fastq2,
        &enriched_scaffolds,
        mapids,
        output_fastq,
        is_paired,
        create_new);
}