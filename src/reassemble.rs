
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::fs::{self,File};
use std::io::{self, BufRead, BufReader, Write};
use std::process::{exit, Command as ProcessCommand, Stdio};
use std::sync::{Arc, RwLock};
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
    bin_qualities: &Arc<RwLock<HashMap<String, BinQuality>>>,
    completeness_cutoff: f64,
    contamination_cutoff: f64
) {
    
    if readfile.is_empty() {
        error!("Error: No read files provided.");
        exit(1);
    }
    if assembler == "spades" {
        // Run SPAdes assembler
        let mut output = ProcessCommand::new("spades.py");
            output.arg("--trusted-contigs")
            .arg(binfile)
            .arg("--only-assembler")
            // .arg("--careful")
            .arg("-o")
            .arg(outputdir)
            .arg("-t")
            .arg(threads.to_string())
            .arg("-m")
            .arg(128.to_string())
            .stdout(Stdio::null());
        if is_paired {
            output.arg("--12").arg(&readfile[0]);
        } else {
            output.arg("-1").arg(&readfile[0]).arg("-2").arg(&readfile[1]);
        }
    
        match output.status() {
            Ok(status) if status.success() => {
                let _ = filterscaffold(&outputdir.join("scaffolds.fasta"));
                let merged_checkm2_output = match assess_bins(
                    &outputdir.join("scaffolds_filtered.fasta"),
                    &outputdir.join("checkm2_results"),
                    threads, 
                    &"fasta"){
                        Ok(path) => path,
                        Err(e) => {
                            error!("Error running assess_bins: {}", e);
                            return;                        }
                    };
                let mergedbin_quality = parse_bins_quality(
                    &merged_checkm2_output
                );
                match mergedbin_quality {
                    Ok(bin_quality_map) => {
                        // Assuming there's only one key-value pair in the HashMap
                        if let Some((_, bin_quality)) = bin_quality_map.into_iter().next() {
                            if bin_quality.contamination < contamination_cutoff && bin_quality.completeness >= completeness_cutoff {
                                let _ = fs::rename(
                                outputdir.join("scaffolds_filtered.fasta"), 
                                resultdir.join(format!("{}_merged.fasta", id)));

                                let mut qualities = bin_qualities
                                    .write()
                                    .unwrap_or_else(|poisoned| poisoned.into_inner());
                                qualities.insert(
                                    format!("{}_merged", id),
                                    BinQuality {
                                        completeness: bin_quality.completeness,
                                        contamination: bin_quality.contamination,
                                    }
                                );

                            } else {
                                let _ = select_highcompletebin(
                                    &component,
                                    bin_qualities,
                                    &bindir,
                                    &resultdir,
                                    completeness_cutoff
                                );
                            }
                        }
                    }
                    Err(e) => {
                        error!("Failed to parse bin qualities: {:?}", e);
                    }
                }
            }
            Ok(_) => {
                warn!("SPAdes failed for {:?} due to low k-mer counts. Selecting best bin."
                    ,binfile.iter()
                    .rev()
                    .nth(2)  // nth(2) gives the third component (0-based index)
                    .map(|bin| bin.to_string_lossy().to_string()).unwrap());
                let _ = select_highcompletebin(
                    &component,
                    bin_qualities,
                    &bindir,
                    &resultdir,
                    completeness_cutoff
                );
            }
            Err(e) => {
                error!("Error: Failed to execute SPAdes command - {}", e);
            }
        }
    } 
    if assembler == "megahit" {
        // Run Megahit assembler
        let mut output = ProcessCommand::new("megahit");
            output.arg("-o")
            .arg(outputdir)
            .arg("-t")
            .arg(threads.to_string())
            .arg("--min-contig-len")
            .arg(500.to_string())
            .stdout(Stdio::null());
        if is_paired {
            output.arg("--12").arg(&readfile[0]);
        } else {
            output.arg("-1").arg(&readfile[0]).arg("-2").arg(&readfile[1]);
        }
    
        match output.status() {
            Ok(status) if status.success() => {
                
                let merged_checkm2_output = match assess_bins(
                    &outputdir.join("final.contigs.fa"),
                    &outputdir.join("checkm2_results"),
                    threads, 
                    &"fa"){
                        Ok(path) => path,
                        Err(e) => {
                            error!("Error running assess_bins: {}", e);
                            return;                        }
                    };
                let mergedbin_quality = parse_bins_quality(
                    &merged_checkm2_output
                );
                match mergedbin_quality {
                    Ok(bin_quality_map) => {
                        // Assuming there's only one key-value pair in the HashMap
                        if let Some((_, bin_quality)) = bin_quality_map.into_iter().next() {
                            if bin_quality.contamination < contamination_cutoff && bin_quality.completeness >= completeness_cutoff {
                                let _ = fs::rename(
                                    outputdir.join("final.contigs.fa"), 
                                    resultdir.join(format!("{}_merged.fasta", id)));
                                let mut qualities = bin_qualities
                                    .write()
                                    .unwrap_or_else(|poisoned| poisoned.into_inner());
                                qualities.insert(
                                    format!("{}_merged", id),
                                    BinQuality {
                                        completeness: bin_quality.completeness,
                                        contamination: bin_quality.contamination,
                                    }
                                );

                            } else {
                                let _ = select_highcompletebin(
                                    &component,
                                    &bin_qualities,
                                    &bindir,
                                    &resultdir,
                                    completeness_cutoff
                                );
                            }
                        }
                    }
                    Err(e) => {
                        error!("Failed to parse bin qualities: {:?}", e);
                    }
                }

            }
            Ok(_) => {
                error!("NOTE: Megahit failed for {:?}! Why?"
                    ,binfile.iter()
                    .rev()
                    .nth(2)  // nth(2) gives the third component (0-based index)
                    .map(|bin| bin.to_string_lossy().to_string()).unwrap());
                let _ = select_highcompletebin(
                    &component,
                    &bin_qualities,
                    &bindir,
                    &resultdir,
                    completeness_cutoff
                );
            }
            Err(e) => {
                error!("Error: Failed to execute MEGAHIT command - {}", e);
            }
        }
    }
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