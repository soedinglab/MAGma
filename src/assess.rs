use csv::ReaderBuilder;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::fs::{self, File};
use std::io::{self};
use std::process::{Command as ProcessCommand, Stdio};
use log::error;

#[derive(Clone)]
pub struct BinQuality {
    pub completeness: f64,
    pub contamination: f64,
}

/// Run CheckM2 to obtain completeness and contamination of input bins
pub fn assess_bins(
    bindir: &PathBuf,
    bincheckm2dir: &PathBuf,
    threads: usize,
    format: &str,
) -> Result<PathBuf, io::Error> {

    let checkm2_qualities = Path::new(bincheckm2dir).join("quality_report.tsv");

    if !checkm2_qualities.exists() {
        let mut output = ProcessCommand::new("checkm2");
        output
        .arg("predict")
        .arg("-i")
        .arg(bindir)
        .arg("-o")
        .arg(bincheckm2dir)
        .arg("-t")
        .arg(threads.to_string())
        .arg("-x")
        .arg(format)
        .arg("--force")
        .stdout(Stdio::null())
        .stderr(Stdio::null());
    
        match output.status() {
            Ok(_) => {
            }
            Err(e) => {
                error!("Error: Failed to execute CheckM2 command - {}. Check if CheckM2 is executable currently", e);
            }
        }
    }    
    Ok(checkm2_qualities)
}

/// Parse CheckM2 result
pub fn parse_bins_quality(
    checkm2_qualities: &PathBuf,
) -> io::Result<HashMap<String, BinQuality>> {

    let _ = File::open(checkm2_qualities).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!(
                "Failed to open checkm2 quality file for bin: {:?}",
                e
            ),
        )
    })?;

    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(File::open(checkm2_qualities)?);
    let mut bin_qualities: HashMap<String, BinQuality> = HashMap::new();

    for result in rdr.records() {
        let record = result?; // Get the record from the CSV
        if record.len() < 3 {
            error!("Skipping invalid record: {:?}", record);
            continue; // Skip records that do not have enough columns
        }
        let bin_id: String = record[0].to_string();
        let completeness: f64 = record[1].parse().unwrap_or(0.0);
        let contamination: f64 = record[2].parse().unwrap_or(0.0);
        bin_qualities.insert(bin_id, 
        BinQuality{ completeness, contamination });
    }
    Ok(bin_qualities)
}

/// Select a high-quality bin if it exist in the cluster
pub fn check_high_quality_bin(
    comp: &HashSet<String>,
    bin_qualities: &HashMap<String, BinQuality>,
    bindir: &PathBuf,
    resultdir: &PathBuf,
    format: &String,
) -> bool {

    let comp_binqualities: HashMap<String, BinQuality> = bin_qualities
        .iter()
        .filter(|(bin, _)| comp.contains(bin.as_str()))  // Compare as &str
        .map(|(bin, q)| (bin.clone(), q.clone()))
        .collect();

    if let Some((binname, _)) = comp_binqualities
        .iter()
        .filter(|(_, q)| q.completeness > 90.0)
        .max_by(|a, b| {
            // select the best bin by quality score
            let quality_score_a = a.1.completeness - (5.0 * a.1.contamination);
            let quality_score_b = b.1.completeness - (5.0 * b.1.contamination);

            quality_score_a
                .partial_cmp(&quality_score_b)
                .unwrap_or(std::cmp::Ordering::Equal) 
                .then_with(|| a.1.contamination.partial_cmp(&b.1.contamination).unwrap()) // Tie-breaker based on contamination
        })

    {
        let bin_path = bindir.join(format!("{}.{}", binname, format));
        let final_path = resultdir.join(format!("{}.fasta", binname));

        if let Err(_) = fs::copy(&bin_path, &final_path) {
            return false;
        }
        return true;
    }
    false
}