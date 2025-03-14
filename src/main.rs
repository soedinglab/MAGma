use std::collections::{HashMap, HashSet};
use std::fs::{self};
use std::io::{self, stderr, Write};
use std::path::PathBuf;
use std::process::exit;
use clap::Parser;
use gfaparser::parse_gfa_fastq;
use rayon::prelude::*;
use rayon::{current_num_threads, ThreadPoolBuilder};
use log::{debug, info};

mod gfaparser;
mod utility;
mod assess;
mod merge;
mod readfetch;
mod reassemble;

// use index::indexfastqreads;
use readfetch::fetch_fastqreads;

fn validate_paths(cli: &Cli) -> io::Result<(PathBuf, Option<PathBuf>, Option<PathBuf>, PathBuf, PathBuf)> {
    let bindir = utility::validate_path(Some(&cli.bindir), "bindir", &cli.format);
    let gfadir = cli.gfadir.as_ref().map(|p| utility::validate_path(
        Some(p), "gfadir", ".gfa").to_path_buf());
    let assemblydir = cli.assemblydir.as_ref().map(|p| utility::validate_path(
        Some(p), "assemblydir", ".fasta").to_path_buf());
    let mapdir = utility::validate_path(Some(&cli.mapdir), "mapdir", "_mapids");
    let readdir = utility::validate_path(Some(&cli.readdir), "readdir", ".fastq");

    Ok((bindir.to_path_buf(), gfadir, assemblydir, mapdir.to_path_buf(), readdir.to_path_buf()))
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Directory containing fasta files of bins
    #[arg(short = 'b', long = "bindir", help = "Directory containing fasta files of bins")]
    bindir: PathBuf,

    /// Average Nucleotide Identity cutoff
    #[arg(short = 'i', long = "ani", default_value_t = 99.0, help = "ANI for clustering bins (%)")]
    ani: f64,

    /// Completeness
    #[arg(short = 'c', long = "completeness", default_value_t = 50.0, help = "Minimum completeness of bins (%)")]
    completeness_cutoff: f64,

    /// Purity
    #[arg(short = 'p', long = "purity", default_value_t = 95.0, help = "Mininum purity of bins (%)")]
    purity_cutoff: f64,

    /// Directory containing mapids files derived from alignment sam/bam files
    #[arg(short = 'm', long = "mapdir", help = "Directory containing mapids files")]
    mapdir: PathBuf,

    /// Directory containing read files
    #[arg(short = 'r', long = "readdir", help = "Directory containing read files")]
    readdir: PathBuf,

    /// Bin file extension
    #[arg(short = 'f', long = "format", default_value = "fasta", help = "Bin file extension")]
    format: String,

    /// Number of threads to use
    #[arg(short = 't', long = "threads", default_value_t = 8, help = "Number of threads to use")]
    threads: usize,

    /// First split bins before merging (if provided, set to true)
    #[arg(long = "split", help = "Split clusters into sample-wise bins before processing")]
    split: bool,
    
    /// Directory containing gfa files for metagenomic samples (in gfa1.2 format)
    #[arg(short = 'g', long = "gfadir", help = "Directory containing gfa files")]
    gfadir: Option<PathBuf>,

    /// Directory containing sample-wise assembly contigs file in fasta format
    #[arg(short = 'a', long = "assemblydir", help = "Directory containing assembly contigs")]
    assemblydir: Option<PathBuf>,

    /// Assembler choice
    #[arg(long = "assembler", default_value = "spades", help = "assembler choice for reassembly step (spades|megahit)")]
    assembler: String,

}

fn main() -> io::Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    // Parse arguments
    let (mut bindir, gfadir, assemblydir, mapdir, readdir) = validate_paths(&cli)?;
    let ani_cutoff = cli.ani;
    let completeness_cutoff = cli.completeness_cutoff;
    let purity_cutoff = cli.purity_cutoff;
    let contamination_cutoff = 100.0 - purity_cutoff;
    let format = cli.format;
    let threads = cli.threads;
    let split = cli.split;
    let assembler: String = cli.assembler;
    let parentdir = bindir.parent().map(PathBuf::from).unwrap_or_else(|| bindir.clone());
    
    println!("Starting merge process with the following parameters:");
    println!("  bindir: {:?}", bindir);
    println!("  ANI (%): {:?}", ani_cutoff);
    println!("  Completeness (%): {:?}", completeness_cutoff);
    println!("  Purity/Contamination (%): {:?}/{:?}", purity_cutoff, contamination_cutoff);
    println!("  gfadir: {:?}", gfadir);
    println!("  assemblydir: {:?}", assemblydir);
    println!("  mapdir: {:?}", mapdir);
    println!("  readdir: {:?}", readdir);
    println!("  format: {}", format);
    println!("  threads: {}", threads);
    println!("  assembler choice: {}", assembler);

    if assembler != "spades" && assembler != "megahit" {
        eprintln!("Error: Invalid assembler choice '{}'. Allowed options: 'spades' or 'megahit'.", assembler);
        exit(1);
    }
    let mut gfa_flag: bool = true;
    let mut gfapath: PathBuf = PathBuf::new();
    if gfadir.is_none() {
        gfa_flag = false;
    } else {
        gfapath = gfadir.unwrap();
    }

    let mut assemblypath: PathBuf = PathBuf::new();
    if !assemblydir.is_none() {
        assemblypath = assemblydir.unwrap();
    }

    let is_paired: bool = utility::check_paired_reads(&readdir);
    if is_paired {
        info!("Detected paired end \
        reads in separate files as \
        <sampleid>_1.fastq[.gz] \
        and <sampleid>_2.fastq[.gz]. \
        They will be processed accordingly.")
    } else {
        info!("Read files are interleaved \
        or single-end as <sampleid>.fastq[.gz]. \
        They will be processed accordingly.")
    }
    
    let binfiles = utility::get_binfiles(&bindir,&format)?;

    if binfiles.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No bin files found. \
            Please provide the correct format argument.",
        ));
    }
   
    // Final merge bins directory
    // eg: resultspath = <bindir>/mergedbins/
    let resultdir: PathBuf = parentdir
            .join("mergedbins");
    if resultdir.exists() {
        fs::remove_dir_all(&resultdir)?;
    }
    fs::create_dir(&resultdir)?;
    
    info!("{} cores are used", current_num_threads());

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
                    debug!("Bin file: {:?}", bin);
                    true
                } else {
                    eprintln!("Bin file does not exist: {:?}", bin);
                    false
                }
            })
            .for_each(|bin| {
                // eg: bin_name = bin_1 (input: <bindir>/bin_1.fa)
                let bin_name = bin.file_stem()
                    .and_then(|name| name.to_str())
                    .unwrap_or_default();
                info!("Splitting bin: {:?}", bin);
                let _ = utility::splitbysampleid(&bin, bin_name, &samplewisebinspath, &format);
            });
        });
        bindir = samplewisebinspath;
        debug!("splitting bins by sample {:?} is completed", bindir);
    }

    // Get sample list
    let bin_sample_map: HashMap<String, String> = utility::get_sample_names(&bindir,&format)?;
    let sample_list: HashSet<String> = bin_sample_map.values().cloned().collect();
    info!("{:?} bin files and {:?} samples found", binfiles.len(), sample_list.len());
    debug!("samples: {:?}", sample_list);

    // Obtain quality of bins
    // eg: checkm2_outputpath = <bindir>/checkm2_results/
    let checkm2_outputpath: PathBuf = bindir
        .join("checkm2_results");
    
    debug!("checkm2 evaluation for {:?} is started", bindir);
    // TODO: remove checkm2 output folder
    let checkm2_qualities = assess::assess_bins(
        &bindir,
        &checkm2_outputpath,
        threads,
        &format)?;

    let bin_qualities = match assess::parse_bins_quality(
        &checkm2_qualities,
    ) {
        Ok(quality) => quality,
        Err(_) => {
            eprintln!(
                "Failed to parse quality inputbins {:?}.",
                bindir
            );
            return Ok(());
        }
    };
    debug!("checkm2 evaluation for {:?} is completed", bindir);

    // None of the bins are pure
    if bin_qualities.len() == 0 {
        debug!("{:?} doesn't have high pure bins. Or if existing checkm2 result is empty, first remove them before running mergebins", bindir);
        return Ok(());
    }        

    let (graph, ani_details) = match merge::calc_ani(&bindir, &bin_qualities, &format, ani_cutoff, contamination_cutoff) {
        Ok((graph, ani_details)) => {
            debug!("Graph was constructed for {:?}", bindir);
            (graph, ani_details)
        },
        Err(e) => {
            eprintln!("Error calculating ANI: {}", e);
            return Ok(()); // Return early in case of an error
        }
    };
    
    debug!("Graph {:?}", graph);
    
    // Assuming `get_connected_samples` or a similar function takes the graph and the ani_details.
    let connected_bins: Vec<HashSet<String>> = merge::get_connected_samples(&graph, &ani_details, ani_cutoff);
    
    pool.install(|| {
        connected_bins
        .par_iter()
        .enumerate()
        .try_for_each(|(id, component)| {
            // Flush stderr once before processing starts
            stderr().flush().ok();

            // Process each connected component
            process_components(
                component.clone(),
                &bindir,
                &gfapath,
                gfa_flag,
                &assemblypath,
                &mapdir,
                &readdir,
                &resultdir,
                bin_sample_map.clone(),
                &format,
                bin_qualities.clone(),
                is_paired,
                assembler.clone(),
                completeness_cutoff,
                contamination_cutoff,
                id,
            )
            .map_err(|e| {
                eprintln!("Error processing bin {:?}: {}", component, e);
                e
            })
        })
        .expect("Error during processing components");
    });

    // clear tmp files
    let _ = fs::remove_dir_all(checkm2_outputpath);
    info!("Bin merging was successfully completed!");  

    Ok(())
}

fn process_components(
    component: HashSet<String>,
    bindir: &PathBuf,
    gfadir: &PathBuf,
    gfa_flag: bool,
    assemblydir: &PathBuf,
    mapdir: &PathBuf,
    readdir: &PathBuf,
    resultdir: &PathBuf,
    bin_sample_map: HashMap<String,String>,
    format: &String,
    bin_qualities: HashMap<String, assess::BinQuality>,
    is_paired: bool,
    assembler:String,
    completeness_cutoff: f64,
    contamination_cutoff: f64,
    id: usize,
) -> std::io::Result<()> {

    // eg: comp = {"binname_S1", "binname_S2"}
    debug!("{} id with length {} for component {:?}", id, component.len(), component);
    if component.len() == 1 {
        debug!("Single component. Checking its quality {:?}", component);
        let binname = component.into_iter().next().expect("The component is empty.");
        if let Some(quality) = bin_qualities.get(&binname) {
            if quality.completeness >= completeness_cutoff {
                let bin_path = bindir.join(format!("{}.{}", binname, format));
                let final_path = resultdir.join(format!("{}.fasta", binname));
                debug!("Copying {:?} to {:?} for component {}", bin_path, final_path, id.to_string());
                
                if let Err(e) = fs::copy(&bin_path, &final_path) {
                    debug!("Failed to copy bin {}: {:?}", binname, e);
                }
            }
        }
        return Ok(());
    }
    
    if assess::check_high_quality_bin(&component, &bin_qualities, bindir, resultdir, id, &format) {
        return Ok(());
    }

    // check quality of the components if merged
    // eg. selected_binset_path = <bindir>/0_combined/
    let selected_binset_path = 
        bindir.join(format!("{}_combined", id.to_string()));
    if selected_binset_path.exists() {
        fs::remove_dir_all(&selected_binset_path)?;
    }
    fs::create_dir(&selected_binset_path)?;

    merge::combine_fastabins(
    &bindir,
    &component,
    &selected_binset_path,
        &format).map_err(|e| {
        eprintln!(
            "Error in combining combined bins for component {}: {}",
            id, e
        );
        e
    })?;

    debug!("input to fetch combined contig ids {}", &selected_binset_path.join("combined.fasta").to_string_lossy());
    let mut all_enriched_scaffolds = HashSet::new();
    if !gfa_flag {
        all_enriched_scaffolds = utility::read_fasta(
            &selected_binset_path.join("combined.fasta").to_string_lossy()
        )?;
    } else {
        let mut create_new = true;
        for samplebin in component.clone() {
            debug!("going to read gfa processing for {}", samplebin);
            let sample = bin_sample_map.get(&samplebin)
                .expect(&format!("Error: File '{}' not found in map!", samplebin));

            // eg: <gfadir>/S1.gfa
            let gfa_path = gfadir.join(format!("{}.gfa", sample));
            let gfa_file = utility::path_to_str(&gfa_path);
            
            // eg: <bindir>/binname_S1.fasta
            let subbin_path = bindir.join(format!("{}.{}", samplebin, format));
            let subbin_file = utility::path_to_str(&subbin_path);
            
            // eg: <assemblydir>/S1.fasta
            let assembly_path = assemblydir.join(format!("{}.fasta", sample));
            let assembly_file = utility::path_to_str(&assembly_path);
            // eg: outputbin = <bindir>/0_combined/S1.fasta
            let enriched_scaffolds = parse_gfa_fastq(
                &sample,
                gfa_file,
                subbin_file,
                assembly_file,
                selected_binset_path
                .join("combined_enriched.fasta"),
                create_new
            )?;
            create_new = false;
            for scaffold in enriched_scaffolds {
                // Check if the scaffold already contains the sample ID separated by 'C'.
                if !scaffold.contains(&format!("{}C", sample)) {
                    // Add sample ID to the scaffold if not present.
                    let modified_scaffold = format!("{}C{}", sample, scaffold);
                    all_enriched_scaffolds.insert(modified_scaffold);
                } else {
                    // If already present, keep the original scaffold.
                    all_enriched_scaffolds.insert(scaffold.clone());
                }
            }
        }
    }
    debug!("Fetched contigs from assembly file");
    debug!("all enriched scaffolds: {} {:?}", all_enriched_scaffolds.len(), all_enriched_scaffolds);
    let scaffold_inputname:&str;
    if gfa_flag {
        scaffold_inputname = "combined_enriched";
    } else {
        scaffold_inputname = "combined";
    }
    for samplebin in component.clone() {
        let sample = bin_sample_map.get(&samplebin)
            .expect(&format!("Error: File '{}' not found in map!", samplebin));


        debug!("processing sample file {} to get reads", sample);
        let mapid_path = mapdir.join(format!("{}_mapids", sample));
        let mapid_file = utility::path_to_str(&mapid_path);
                    
        let read_files: Vec<String> = if is_paired {
            let read_path1 = utility::find_file_with_extension(&readdir, &format!("{}_1", sample));
            let read_path2 = utility::find_file_with_extension(&readdir, &format!("{}_2", sample));
    
            vec![
                read_path1.to_str().expect("Failed to convert PathBuf to &str").to_string(),
                read_path2.to_str().expect("Failed to convert PathBuf to &str").to_string()
            ]
        } else {
            let read_path = utility::find_file_with_extension(&readdir, &sample);
            
            vec![
                read_path.to_str().expect("Failed to convert PathBuf to &str").to_string()
            ]
        };
        
        debug!("map id {} read file1 {:?} ", mapid_file, read_files);
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
        // For paired-end reads, return both "_1" and "_2" fastq files
        vec![
            selected_binset_path.join(format!("{}_1.fastq", scaffold_inputname)),
            selected_binset_path.join(format!("{}_2.fastq", scaffold_inputname)),
        ]
    } else {
        // For single-end reads, return just the "combined.fastq" file
        vec![selected_binset_path.join(format!("{}.fastq", scaffold_inputname))]
    };

    info!("Performing assembly for {} with reads {:?}", id, reads_path);
    let reassembly_outputdir = selected_binset_path.join("assembly");

    let _ = reassemble::run_reassembly(
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
        completeness_cutoff,
        contamination_cutoff,
    );
    info!("Component {} is processed", id);
    
    // clean folders
    let _ = fs::remove_dir_all(selected_binset_path);
    Ok(())
}
