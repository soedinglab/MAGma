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

pub fn assess_bins(
    bindir: &PathBuf,
    bincheckm2: &PathBuf,
    threads: usize,
    format: &str,
) -> Result<PathBuf, io::Error> {
    let checkm2_qualities = Path::new(bincheckm2).join("quality_report.tsv");
    // Run CheckM2 run
    if !checkm2_qualities.exists() {
        // println!("{:?}/quality_report.tsv not found. Running checkm2...", bindir);
        let mut output = ProcessCommand::new("checkm2");
        output
        .arg("predict")
        .arg("-i")
        .arg(bindir)
        .arg("-o")
        .arg(bincheckm2)
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
    // read checkm2 output file
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

pub fn select_highcompletebin(
    bin_samplenames: &HashSet<String>,
    bin_qualities: &HashMap<String, BinQuality>,
    bindir: &PathBuf,
    outputpath: &PathBuf,
    completeness_cutoff: f64,
) -> io::Result<()> {
    
    let highest_completebin = bin_samplenames
        .iter()
        .filter_map(|bin| bin_qualities.get(bin).map(|quality| (bin, quality.completeness)))
        .max_by(|(_, completeness1), (_, completeness2)| {
            completeness1.partial_cmp(completeness2).unwrap_or(std::cmp::Ordering::Equal)
        })
        .filter(|(_, max_completeness)| *max_completeness >= completeness_cutoff)
        .map(|(bin, _)| bin.clone()); 

    if let Some(sample_id) = highest_completebin {
        let bin_path = bindir.join(format!("{}.fasta", sample_id));
        
        fs::copy(bin_path, outputpath.join(format!("{}.fasta", sample_id))).ok();
    }
    Ok(())

}

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
        .max_by(|a, b| {a.1.completeness
            .partial_cmp(&b.1.completeness)
            .unwrap()
            .then_with(|| b.1.contamination.partial_cmp(&a.1.contamination).unwrap())
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