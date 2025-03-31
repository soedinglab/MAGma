
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::process::{exit, Command as ProcessCommand, Stdio};
use std::sync::{Arc, Mutex};
use log::{error, warn};
use crate::assess::{assess_bins, parse_bins_quality, BinQuality};

pub fn run_reassembly(
    readfile: &[PathBuf],
    binfile: &PathBuf,
    outputdir: &PathBuf,
    is_paired: bool,
    threads: usize,
    assembler: &String,
    resultdir: &PathBuf,
    id: usize,
    bindir: &PathBuf,
    component: &HashSet<String>,
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

    let mut selected_bin: Option<String> = None;
    let mut selected_completeness: Option<f64> = None;
    let mut selected_contamination: Option<f64> = None;

    if let Some((bin_name, completeness, contamination)) = component
        .iter()
        .filter_map(|bin| {
            bin_qualities.get(bin).map(|quality| (bin, quality.completeness, quality.contamination))
        })
        .filter(|(_, completeness, _)| *completeness >= completeness_cutoff)
        .max_by(|(_, completeness1, contamination1), (_, completeness2, contamination2)| {
            // select the best bin by quality score
            let score1 = completeness1 - (5.0 * contamination1);
            let score2 = completeness2 - (5.0 * contamination2);
            
            score1
                .partial_cmp(&score2)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| contamination1.partial_cmp(contamination2).unwrap_or(std::cmp::Ordering::Equal).reverse())
        })
    {
        selected_bin = Some(bin_name.to_string());
        selected_completeness = Some(completeness);
        selected_contamination = Some(contamination);
    }

    let selected_quality_score = selected_completeness
        .zip(selected_contamination)
        .map(|(completeness, contamination)| completeness - (5.0 * contamination))
        .unwrap_or(completeness_cutoff - (5.0 * contamination_cutoff));

    if let Err(e) = command_status {
        error!("Assembler failed: {}", e);
        let _ = select_highcompletebin(
            selected_bin,
            bindir,
            resultdir,
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
                    let quality_score = bin_quality.completeness - (5.0 * bin_quality.contamination);
                    if bin_quality.contamination < contamination_cutoff
                        && bin_quality.completeness >= completeness_cutoff
                        && quality_score > selected_quality_score

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
                            selected_bin,
                            bindir,
                            resultdir,
                        );
                    }
                }
        } else {
            warn!("Failed to parse bin qualities.");
            let _ = select_highcompletebin(
                selected_bin,
                bindir,
                resultdir,
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
        .stdout(Stdio::null())
        .stderr(Stdio::null());

    if is_paired {
        cmd.arg("-1").arg(&readfile[0]).arg("-2").arg(&readfile[1]);
    } else {
        cmd.arg("-s").arg(&readfile[0]);
    }

    cmd.status().map(|status| {
        if !status.success() {
            Err(std::io::Error::new(std::io::ErrorKind::Other, "SPAdes failed due to low k-mer counts"))
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
        .stdout(Stdio::null())
        .stderr(Stdio::null());

    if is_paired {
        cmd.arg("-1").arg(&readfile[0]).arg("-2").arg(&readfile[1]);
    } else {
        cmd.arg("-r").arg(&readfile[0]);
    }

    cmd.status().map(|status| {
        if !status.success() {
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
            current_sequence.push_str(line.trim());
        }
    }

    if !current_header.is_empty() && current_sequence.len() >= 500 {
        writeln!(output, "{}", current_header)?;
        writeln!(output, "{}", current_sequence)?;
    }
    Ok(())
}

fn select_highcompletebin(
    selected_bin: Option<String>,
    bindir: &PathBuf,
    outputpath: &PathBuf,
) -> io::Result<()> {

    if let Some(sample_id) = selected_bin {
        let bin_path = bindir.join(format!("{}.fasta", sample_id));
        fs::copy(bin_path, outputpath.join(format!("{}.fasta", sample_id))).ok();
    }
    Ok(())

}