use std::collections::{HashMap, HashSet};
use std::fs::{self};
use std::io::{self, stderr, Write};
use std::path::PathBuf;
use std::process::exit;
use assess::BinQuality;
use clap::Parser;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use log::{debug, error, info, warn};
use std::sync::{Arc, Mutex};
use readfetch::fetch_fastqreads;

mod utility;
mod assess;
mod merge;
mod readfetch;
mod reassemble;

// check for valid input paths
fn validate_paths(cli: &Cli) -> io::Result<(PathBuf, PathBuf, PathBuf)> {
    let bindir = utility::validate_path(Some(&cli.bindir), "bindir", &cli.format);
    let mapdir = utility::validate_path(Some(&cli.mapdir), "mapdir", "_mapids");
    let readdir = utility::validate_path(Some(&cli.readdir), "readdir", ".fastq");

    Ok((bindir.to_path_buf(), mapdir.to_path_buf(), readdir.to_path_buf()))
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Directory containing fasta files of bins
    #[arg(short = 'b', long = "bindir", help = "Directory containing fasta files of bins")]
    bindir: PathBuf,

    /// Directory containing read files
    #[arg(short = 'r', long = "readdir", help = "Directory containing read files")]
    readdir: PathBuf,

    /// Directory containing mapids files derived from alignment sam/bam files
    #[arg(short = 'm', long = "mapdir", help = "Directory containing mapids files")]
    mapdir: PathBuf,

    /// Average Nucleotide Identity cutoff
    #[arg(short = 'i', long = "ani", default_value_t = 99.0, help = "ANI for clustering bins (%)")]
    ani: f64,

    /// Completeness
    #[arg(short = 'c', long = "completeness", default_value_t = 50.0, help = "Minimum completeness of bins (%)")]
    completeness_cutoff: f64,

    /// Purity
    #[arg(short = 'p', long = "purity", default_value_t = 95.0, help = "Mininum purity (1- contamination) of bins (%)")]
    purity_cutoff: f64,

    /// Bin file extension
    #[arg(short = 'f', long = "format", default_value = "fasta", help = "Bin file extension")]
    format: String,

    /// Number of threads to use
    #[arg(short = 't', long = "threads", default_value_t = 8, help = "Number of threads to use")]
    threads: usize,

    /// First split bins before merging (if provided, set to true)
    #[arg(long = "split", help = "Split clusters into sample-wise bins before processing")]
    split: bool,

    /// CheckM2 quality file
    #[arg(short = 'q', long = "qual", help = "Quality file produced by CheckM2 (quality_report.tsv)")]
    qual: Option<PathBuf>,

    /// Assembler choice
    #[arg(long = "assembler", default_value = "spades", help = "assembler choice for reassembly step (spades|megahit), spades is recommended")]
    assembler: String,

}

fn main() -> io::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let cli = Cli::parse();

    // Parse arguments
    let (mut bindir, mapdir, readdir) = validate_paths(&cli)?;
    let ani_cutoff = cli.ani;
    let completeness_cutoff = cli.completeness_cutoff;
    let purity_cutoff = cli.purity_cutoff;
    let contamination_cutoff = 100.0 - purity_cutoff;
    let format = cli.format;
    let threads = cli.threads;
    let split = cli.split;
    let assembler: String = cli.assembler;
    let qual = cli.qual;
    let parentdir = bindir.parent().map(PathBuf::from).unwrap_or_else(|| bindir.clone());
    
    info!("Starting MAGma with parameters:");
    info!("  🔹 Bins Directory: {:?}", bindir);
    info!("  🔹 ANI Cutoff: {:.1}%", cli.ani);
    info!("  🔹 Completeness Cutoff: {:.1}%", cli.completeness_cutoff);
    info!("  🔹 Purity/Contamination: {:.1}%/{:.1}%", cli.purity_cutoff, contamination_cutoff);
    info!("  🔹 Map Directory: {:?}", mapdir);
    info!("  🔹 Read Directory: {:?}", readdir);
    info!("  🔹 File Format: {}", format);
    info!("  🔹 Threads: {}", threads);
    info!("  🔹 Assembler: {}", assembler);

    if !["spades", "megahit"].contains(&assembler.as_str()) {
        error!("Error: Invalid assembler choice '{}'. Allowed options: 'spades' or 'megahit'.", assembler);
        exit(1);
    }

    let is_paired: bool = utility::check_paired_reads(&readdir);
    if is_paired {
        info!("Detected paired end \
        reads in separate files as \
        <sampleid>_1.fastq \
        and <sampleid>_2.fastq.")
    } else {
        info!("Detected single-end reads as <sampleid>.fastq.")
    }
    
    let binfiles = utility::get_binfiles(&bindir,&format)?;

    if binfiles.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No bin files found. \
            Please provide the correct format argument.",
        ));
    }
   
    // Output directory
    // eg: resultspath = <parentpathof_bindir>/mags_90comp_95purity/
    let resultdir: PathBuf = parentdir
        .join(format!(
        "mags_{}comp_{}purity",
        completeness_cutoff as u32,
        purity_cutoff as u32
    ));
    if resultdir.exists() {
        info!("Output folder: {:?} already exist. Cleaning it", &resultdir);
        fs::remove_dir_all(&resultdir)?;
    }
    fs::create_dir(&resultdir)?;
    
    let pool = ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Failed to build thread pool");

    // Split bin by sample id
    if split {
        // Create a directory to store sample-wise bins
        let samplewisebinspath: PathBuf = parentdir
            .join("samplewisebins");
        if samplewisebinspath.exists() {
            fs::remove_dir_all(&samplewisebinspath)?;
        }
        fs::create_dir(&samplewisebinspath)?;

        pool.install(|| {
        binfiles.par_iter()
            .filter_map(|bin| bin.canonicalize().ok())
            .filter(|bin| {
                if bin.exists() {
                    true
                } else {
                    error!("Bin file does not exist: {:?}", bin);
                    false
                }
            })
            .for_each(|bin| {
                // eg: bin_name = bin_1 (input: <bindir>/bin_1.fa)
                let bin_name = bin.file_stem()
                    .and_then(|name| name.to_str())
                    .unwrap_or_default();
                utility::splitbysampleid(&bin, bin_name, &samplewisebinspath, &format).ok();
            });
        });
        bindir = samplewisebinspath;
        info!("splitting bins by sample {:?} is completed", bindir);
    }

    // Get sample list
    let bin_sample_map: HashMap<String, String> = utility::get_sample_names(&bindir,&format)?;
    let sample_count= bin_sample_map.values().collect::<HashSet<_>>().len();
    info!("{:?} bin files and {:?} samples found", binfiles.len(), sample_count);

    // Obtain quality of bins
    // eg: checkm2_outputpath = <parentpathof_bindir>/mags_90comp_95purity/checkm2_results/
    let checkm2_outputpath: PathBuf = resultdir
        .join("checkm2_results");
    debug!("checkm2 output path {:?}",checkm2_outputpath);
    let checkm2_qualities = if let Some(qual_path) = &qual {
        // User have alredy provided CheckM2 quality file
        if qual_path.is_file() && fs::metadata(qual_path).map(|m| m.len() > 0).unwrap_or(false) {
            qual_path.clone()
        } else {
            info!("Provided quality file {:?} is missing or empty. Running CheckM2...", qual_path);
            assess::assess_bins(
                &bindir,
                &checkm2_outputpath,
                threads,
                &format
            )
            .expect("Failed to run CheckM2")
        }
    } else {
        assess::assess_bins(
            &bindir,
            &checkm2_outputpath,
            threads,
            &format
        )
        .expect("Failed to run CheckM2")
    };

    // Obtain bins quality and store in a hashmap
    let mut bin_qualities = match assess::parse_bins_quality(
        &checkm2_qualities,
    ) {
        Ok(quality) => quality,
        Err(_) => {
            error!(
                "Failed to parse CheckM2 quality of inputbins {:?}. Check input --format option and if DIAMOND database is accessible for CheckM2",
                bindir
            );
            return Ok(());
        }
    };
    info!("checkm2 evaluation for {:?} is completed", bindir);

    // None of the bins are pure
    if bin_qualities.len() == 0 {
        info!("Input {:?} does not have any high pure bins. Or if existing checkm2 result is empty, first remove them before running maga", bindir);
        return Ok(());
    }        

    debug!("Bin qualities length before reassembly: {}", bin_qualities.len());
    let (graph, ani_details) = match merge::calc_ani(&bindir, &bin_qualities, &format, ani_cutoff, contamination_cutoff) {
        Ok((graph, ani_details)) => {
            (graph, ani_details)
        },
        Err(e) => {
            error!("Error calculating ANI: {}", e);
            return Ok(());
        }
    };
    
    // Cluster bins based on ANI
    let connected_bins: Vec<HashSet<String>> = merge::get_connected_samples(&graph, &ani_details, ani_cutoff);
    
    // Collect completeness and purity of merged and reassembled bins
    let merged_bin_qualities: Arc<Mutex<HashMap<String, BinQuality>>> = Arc::new(Mutex::new(HashMap::new()));

    pool.install(|| {
        connected_bins
        .par_iter()
        .enumerate()
        .try_for_each(|(id, component)| {
            // Flush stderr once before processing starts
            stderr().flush().ok();
            // Process each connected component
            let merged_bin_quality = Arc::clone(&merged_bin_qualities);
            process_components(
                &component,
                &bindir,
                &mapdir,
                &readdir,
                &resultdir,
                &bin_sample_map,
                &format,
                &bin_qualities,
                &merged_bin_quality,
                is_paired,
                &assembler,
                completeness_cutoff,
                contamination_cutoff,
                id,
            )
            .map_err(|e| {
                error!("Error processing bin {:?}: {}", component, e);
                e
            })
        })
        .expect("Error during processing components");
    });
    
    {
        let merged_bin_qualities = merged_bin_qualities.lock().unwrap();

        for (key, value) in merged_bin_qualities.iter() {
            bin_qualities.insert(key.to_string(), value.clone());
        }
    }

    // Final dereplication using skani
    let _ = merge::drep_finalbins(&resultdir, &bin_qualities, ani_cutoff);
       
    info!("MAGma is successfully completed!");  

    Ok(())
}

/// Process cluster in parallel
fn process_components(
    component: &HashSet<String>,
    bindir: &PathBuf,
    mapdir: &PathBuf,
    readdir: &PathBuf,
    resultdir: &PathBuf,
    bin_sample_map: &HashMap<String,String>,
    format: &String,
    bin_qualities: &HashMap<String, BinQuality>,
    merged_bin_quality: &Arc<Mutex<HashMap<String, BinQuality>>>,
    is_paired: bool,
    assembler: &String,
    completeness_cutoff: f64,
    contamination_cutoff: f64,
    id: usize,
) -> io::Result<()> {

    debug!("{} id component with bins {:?}", id, component);
    
    // eg: comp = {"binname_S1", "binname_S2"}

    // Singleton cluster, save the bin in the output
    if component.len() == 1 {
        let binname = component.into_iter().next().expect("The component is empty.");    

        if let Some(quality) = bin_qualities.get(binname) {
            if quality.completeness >= completeness_cutoff {
                let bin_path = bindir.join(format!("{}.{}", binname, format));
                let final_path = resultdir.join(format!("{}.fasta", binname));
                fs::copy(&bin_path, &final_path).ok();
            }
        }
        return Ok(());
    }
    
    // Check if already cluster has the high-quality bin (>90% comp, <5% cont)
    if assess::check_high_quality_bin(&component, &bin_qualities, bindir, resultdir, &format) {
        return Ok(());
    }

    // eg. selected_binset_path = <bindir>/0_combined/
    let selected_binset_path = 
        resultdir.join(format!("{}_combined", id));
    if selected_binset_path.exists() {
        fs::remove_dir_all(&selected_binset_path)?;
    }
    fs::create_dir(&selected_binset_path)?;

    // Merge bins within the cluster
    merge::combine_fastabins(
    &bindir,
    &component,
    &selected_binset_path,
        format).map_err(|e| {
        error!(
            "Error in combining combined bins for component {}: {}",
            id, e
        );
        e
    })?;

    // (Obsolete) Enrich bins by adding contigs that are directly linked to bin in the assembly graph
    let all_enriched_scaffolds = utility::read_fasta(
        &selected_binset_path.join("combined.fasta").to_string_lossy()
    )?;
    

    let scaffold_inputname:&str = "combined";


    // Collect reads mapped to contigs in the merged set
    for samplebin in component {
        let sample = bin_sample_map.get(samplebin)
            .unwrap_or_else(|| panic!("Error: File '{}' not found in map!", samplebin));

        let mapid_path = mapdir.join(format!("{}_mapids", sample));
        let mapid_file = utility::path_to_str(&mapid_path);
                    
        let read_files: Vec<String> = if is_paired {
            let read_path1 = utility::find_file_with_extension(readdir, &format!("{}_1", sample));
            let read_path2 = utility::find_file_with_extension(readdir, &format!("{}_2", sample));
    
            vec![
                read_path1.to_str().expect("Failed to convert PathBuf to &str").to_string(),
                read_path2.to_str().expect("Failed to convert PathBuf to &str").to_string()
            ]
        } else {
            let read_path = utility::find_file_with_extension(readdir, sample);
            
            vec![
                read_path.to_str().expect("Failed to convert PathBuf to &str").to_string()
            ]
        };
        
        let _ = fetch_fastqreads(
            &all_enriched_scaffolds,
            mapid_file,
            read_files,
            selected_binset_path
            .join(format!("{}.fasta", scaffold_inputname)),
            is_paired,
        );
    }
    let reads_path = if is_paired {
        // For paired-end reads, return both "combined_1" and "combined_2" fastq files
        vec![
            selected_binset_path.join(format!("{}_1.fastq", scaffold_inputname)),
            selected_binset_path.join(format!("{}_2.fastq", scaffold_inputname)),
        ]
    } else {
        // For single-end reads, return just the "combined.fastq" file
        vec![selected_binset_path.join(format!("{}.fastq", scaffold_inputname))]
    };

    let reassembly_outputdir = selected_binset_path.join("assembly");

    // Reassemble merged bins using contigs and mapped reads by SPAdes
    reassemble::run_reassembly(
        &reads_path,
        &selected_binset_path.join(format!("{}.fasta", scaffold_inputname)),
        &reassembly_outputdir,
        true,
        8,
        assembler,
        resultdir,
        id,
        bindir,
        component,
        bin_qualities,
        merged_bin_quality,
        completeness_cutoff,
        contamination_cutoff,
        format
    );
    info!("Reassembly is completed for component {}", id.to_string());
    
    // clean folder
    if !cfg!(debug_assertions) {
        if let Err(e) = fs::remove_dir_all(&selected_binset_path) {
            warn!("Failed to delete folder {:?}: {}", selected_binset_path, e);
        }
    }
    Ok(())
}
