
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::fs::{self,File};
use std::io::{self, BufRead, BufReader, Write};
use std::process::{exit, Command as ProcessCommand, Stdio};
use std::sync::{Arc, Mutex};
use log::{error, warn};
use crate::assess::{assess_bins, parse_bins_quality, select_highcompletebin, BinQuality};

pub fn run_reassembly(
    readfile: &[PathBuf],
    binfile: &PathBuf,
    outputdir: &PathBuf,
    is_paired: bool,
    threads: usize,
    assembler: String,
    resultdir: &PathBuf,
    id: usize,
    bindir: &PathBuf,
    component: HashSet<String>,
    bin_qualities: &HashMap<String, BinQuality>,
    merged_bin_quality: &Arc<Mutex<HashMap<String, BinQuality>>>,
    completeness_cutoff: f64,
    contamination_cutoff: f64
) {
    
    if readfile.is_empty() {
        error!("Error: No read files provided.");
        exit(1);
    }

    let command_status = match assembler.as_str() {
        "spades" => {
            let status = run_spades(readfile, binfile, outputdir, is_paired, threads);
            if status.is_ok() {
                let _ = filterscaffold(&outputdir.join("scaffolds.fasta")); // Added filtering step
            }
            status
        }
        "megahit" => run_megahit(readfile, outputdir, is_paired, threads),
        _ => {
            error!("Unknown assembler: {}", assembler);
            return;
        }
    };

    if let Err(e) = command_status {
        error!("Assembler failed: {}", e);
        let _ = select_highcompletebin(
            &component,
            bin_qualities,
            bindir,
            resultdir,
            completeness_cutoff
        );
        return;
    }

    let merged_bin_path = if assembler == "spades" {
        outputdir.join("scaffolds_filtered.fasta")
    } else {
        outputdir.join("final.contigs.fa")
    };

    if let Ok(merged_checkm2_output) =
        assess_bins(&merged_bin_path, &outputdir.join("checkm2_results"), threads, "fasta")
        {
            if let Ok(mut bin_quality_map) = parse_bins_quality(&merged_checkm2_output) {
                if let Some((_, bin_quality)) = bin_quality_map.drain().next() {
                    if bin_quality.contamination < contamination_cutoff
                        && bin_quality.completeness >= completeness_cutoff
                    {
                        let _ = fs::rename(&merged_bin_path, resultdir.join(format!("{}_merged.fasta", id)));
                        match merged_bin_quality.lock() {
                            Ok(mut mergedbin_quality_map) => {
                                mergedbin_quality_map.insert(format!("{}_merged", id), bin_quality);
                            },
                            Err(e) => {
                                warn!("Error locking the mutex: {}", e);
                            }
                        }
                    } else {
                        let _ = select_highcompletebin(
                            &component,
                            bin_qualities,
                            bindir,
                            resultdir,
                            completeness_cutoff
                        );
                    }
                }
        } else {
            warn!("Failed to parse bin qualities.");
            let _ = select_highcompletebin(
                &component,
                bin_qualities,
                bindir,
                resultdir,
                completeness_cutoff
            );
        }
    }

}

fn run_spades(
    readfile: &[PathBuf],
    binfile: &PathBuf,
    outputdir: &PathBuf,
    is_paired: bool,
    threads: usize,
) -> std::io::Result<()> {
    let mut cmd = ProcessCommand::new("spades.py");
    cmd.arg("--trusted-contigs")
        .arg(binfile)
        .arg("--only-assembler")
        .arg("-o")
        .arg(outputdir)
        .arg("-t")
        .arg(threads.to_string())
        .arg("-m")
        .arg(128.to_string())
        .stdout(Stdio::null());

    if is_paired {
        cmd.arg("--12").arg(&readfile[0]);
    } else {
        cmd.arg("-1").arg(&readfile[0]).arg("-2").arg(&readfile[1]);
    }

    cmd.status().map(|status| {
        if !status.success() {
            error!("SPAdes failed due to low k-mer counts.");
            Err(std::io::Error::new(std::io::ErrorKind::Other, "SPAdes failed"))
        } else {
            Ok(())
        }
    })?
}


fn run_megahit(
    readfile: &[PathBuf],
    outputdir: &PathBuf,
    is_paired: bool,
    threads: usize,
) -> std::io::Result<()> {
    let mut cmd = ProcessCommand::new("megahit");
    cmd.arg("-o")
        .arg(outputdir)
        .arg("-t")
        .arg(threads.to_string())
        .arg("--min-contig-len")
        .arg("500")
        .stdout(Stdio::null());

    if is_paired {
        cmd.arg("--12").arg(&readfile[0]);
    } else {
        cmd.arg("-1").arg(&readfile[0]).arg("-2").arg(&readfile[1]);
    }

    cmd.status().map(|status| {
        if !status.success() {
            error!("MEGAHIT failed.");
            Err(std::io::Error::new(std::io::ErrorKind::Other, "MEGAHIT failed"))
        } else {
            Ok(())
        }
    })?
}

fn filterscaffold(input_file: &PathBuf) -> io::Result<()> {

    let output_file = match Path::new(input_file)
        .extension() {
        Some(ext) => {
            let stem = Path::new(input_file)
                .file_stem()
                .unwrap_or_default();
            let parent = Path::new(input_file)
                .parent()
                .unwrap_or_else(||
                Path::new(""));
            parent.join(
        format!("{}_filtered.{}",
            stem.to_string_lossy(),
            ext.to_string_lossy()))
        }
        None => Path::new(input_file)
            .with_file_name(
            format!("{}_filtered",
            input_file.display())),
    };
    let input = File::open(input_file)?;
    let mut output = File::create(&output_file)?;
    let reader = BufReader::new(input);

    let mut current_header = String::new();
    let mut current_sequence = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // Process the previous sequence, if any.
            if !current_header.is_empty() && current_sequence.len() >= 500 {
                writeln!(output, "{}", current_header)?;
                writeln!(output, "{}", current_sequence)?;
            }
            current_header = line;
            current_sequence.clear();
        } else {
            current_sequence.push_str(&line.trim());
        }
    }

    if !current_header.is_empty() && current_sequence.len() >= 500 {
        writeln!(output, "{}", current_header)?;
        writeln!(output, "{}", current_sequence)?;
    }
    Ok(())
}