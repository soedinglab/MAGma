use std::collections::HashSet;
use std::process::{Command, Stdio};
use std::path::PathBuf;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, Write};
use crate::utility;
use log::error;

pub fn fetch_fastqreads(
    enriched_scaffolds: &HashSet<String>,
    mapids: &str,
    fastq_files: Vec<String>,
    outputbin: PathBuf,
    is_paired: bool
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

    let _ = write_selected_reads(
        fastq_files,
        &enriched_scaffolds,
        mapids,
        &output_fastq,
        is_paired);

    let file_metadata = std::fs::metadata(&output_fastq[0]).map_err(|e| {
        error!("Error accessing {:?}: {}", output_fastq, e);
        Box::new(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("Failed to access {:?}", output_fastq),
        )) as Box<dyn std::error::Error>
    })?;

    if file_metadata.len() == 0 {
        error!("Error: The output file {:?} is empty", output_fastq[0]);
        return Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("No reads were written to {:?}", output_fastq[0]),
        )));
    }
    Ok(())
}

fn write_selected_reads(
    fastq_files: Vec<String>,
    enriched_scaffolds: &HashSet<String>,
    mapid_file: &str,
    output_fastq: &Vec<String>,
    is_paired: bool,
) -> Result<(), io::Error> {
    let mfile = File::open(&mapid_file)?;
    let mapid_reader = io::BufReader::new(mfile);
    let readid_file = format!("{}_readids", output_fastq[0].replace(".fastq", ""));
    let mut idfile = io::BufWriter::new(
        File::create(&readid_file)?);
    for line in mapid_reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            error!("map id file is not in two column format");
            continue;
        }

        // As of now it only works for readid format: @SRR3961047.1
        let read_id = parts[0];
        let scaffold_id = parts[1];

        // Step 4: Check if scaffold_id exists in the enriched set
        if enriched_scaffolds.contains(scaffold_id) {
            writeln!(idfile, "{}", read_id)?;
        }
        idfile.flush()?;
    }

    let mut outfile1 = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&output_fastq[0])?;

    let mut outfile2 = if is_paired {
        Some(OpenOptions::new()
            .create(true)
            .append(true)
            .open(&output_fastq[1])?)
    } else {
        None
    };

    let process_seqtk = |fastq: &String, outfile: &mut File| -> Result<(), io::Error> {
        let mut child = Command::new("seqtk")
            .arg("subseq")
            .arg(fastq)
            .arg(&readid_file)
            .stdout(Stdio::piped())
            .spawn()?;
    
        io::copy(&mut child.stdout.take().unwrap(), outfile)?;
        Ok(())
    };

    process_seqtk(&fastq_files[0], &mut outfile1)?;
    // Process reverse reads
    if let Some(outfile2) = &mut outfile2 {
        process_seqtk(&fastq_files[1], outfile2)?;
    }

    Ok(())
}
