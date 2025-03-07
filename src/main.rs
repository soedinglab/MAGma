use std::collections::{HashMap, HashSet};
use std::fs::{self};
use std::io::{self, stderr, Write};
use std::path::PathBuf;
use std::process::exit;
use clap::Parser;
use gfaparser::parse_gfa_fastq;
use rayon::prelude::*;
use rayon::{current_num_threads, ThreadPoolBuilder};
use log::{debug, error, info};
use num_cpus;

mod gfaparser;
mod readfetch;
mod utility;

// use index::indexfastqreads;
use readfetch::fetch_fastqreads;

fn validate_paths(cli: &Cli) -> io::Result<(PathBuf, Option<PathBuf>, PathBuf, PathBuf, PathBuf)> {
    let bindir = utility::validate_path(Some(&cli.bindir), "bindir", &cli.format);
    let gfadir = cli.gfadir.as_ref().map(|p| utility::validate_path(Some(p), "gfadir", ".gfa").to_path_buf());
    let assemblydir = utility::validate_path(Some(&cli.assemblydir), "assemblydir", ".fasta");
    let mapdir = utility::validate_path(Some(&cli.mapdir), "mapdir", "_mapids");
    let readdir = utility::validate_path(Some(&cli.readdir), "readdir", ".fastq");

    Ok((bindir.to_path_buf(), gfadir, assemblydir.to_path_buf(), mapdir.to_path_buf(), readdir.to_path_buf()))
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Directory containing fasta files of bins
    #[arg(short = 'b', long = "bindir", help = "Directory containing fasta files of bins")]
    bindir: PathBuf,

    /// Average Nucleotide Identity cutoff
    #[arg(short = 'i', long = "ani", default_value_t = 99.0, help = "ANI for clustering bins")]
    ani: f32,

    /// Directory containing gfa files for metagenomic samples (in gfa1.2 format)
    #[arg(short = 'g', long = "gfadir", help = "Directory containing gfa files")]
    gfadir: Option<PathBuf>,

    /// Directory containing sample-wise assembly contigs file in fasta format
    #[arg(short = 'a', long = "assemblydir", help = "Directory containing assembly contigs")]
    assemblydir: PathBuf,

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

    /// Minimum overlap length
    #[arg(short = 'l', long = "min_overlaplen", default_value_t = 1000, help = "Minimum overlap length")]
    min_overlaplen: usize,

    /// First split bins before merging (if provided, set to true)
    #[arg(long = "split", help = "Split clusters into sample-wise bins before processing")]
    split: bool,

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
    let format = cli.format;
    let threads = cli.threads;
    let min_overlaplen = cli.min_overlaplen;
    let split = cli.split;
    let assembler: String = cli.assembler;
    let parentdir = bindir.parent().map(PathBuf::from).unwrap_or_else(|| bindir.clone());
    
    println!("Starting merge process with the following parameters:");
    println!("  bindir: {:?}", bindir);
    println!("  ANI (%): {:?}", ani_cutoff);
    println!("  gfadir: {:?}", gfadir);
    println!("  assemblydir: {:?}", assemblydir);
    println!("  mapdir: {:?}", mapdir);
    println!("  readdir: {:?}", readdir);
    println!("  format: {}", format);
    println!("  threads: {}", threads);
    println!("  min_overlaplen: {}", min_overlaplen);
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

    let sample_list = utility::get_gfa_file_names(&assemblydir);
    info!("{:?} bin files and {:?} samples found", binfiles.len(), sample_list.len());
    
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
    // Obtain quality of bins
    // eg: checkm2_outputpath = <bindir>/checkm2_results/
    let checkm2_outputpath: PathBuf = bindir
        .join("checkm2_results");
    
    debug!("checkm2 evaluation for {:?} is started", bindir);
    // TODO: remove checkm2 output folder
    let checkm2_qualities = utility::assess_bins(
        &bindir,
        &checkm2_outputpath,
        num_cpus::get()-2,
        &format)?;

    let bin_qualities = match utility::parse_bins_quality(
        &checkm2_qualities,
    ) {
        Ok(quality) => quality,
        Err(_) => {
            eprintln!(
                "Failed to parse quality inputbins {:?}.",
                bindir
            );
            return Ok(()); // Skip further processing for this bin
        }
    };
    debug!("checkm2 evaluation for {:?} is completed", bindir);

    // None of the bins are pure
    if bin_qualities.len() == 0 {
        debug!("{:?} doesn't have high pure bins. Or if existing checkm2 result is empty, first remove them before running mergebins", bindir);
        return Ok(());
    }        

    let result = utility::calc_ani(&bindir, &bin_qualities, &format, ani_cutoff);
    debug!("Graph was constructed for {:?}", bindir);
    let connected_bins: Vec<HashSet<String>> = match result {
        Ok(graph) => {
            debug!("Graph {:?}", graph);
            let bins = utility::get_connected_samples(&graph);
            // let bins = utility::find_cliques(&graph);
            debug!("Connected components {:?} for {:?}", bins, bindir);
            bins
        },
        Err(e) => {
            eprintln!("Error calculating ANI: {}", e);
            return Ok(());
        }
    };

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
                &assemblydir,
                &mapdir,
                &readdir,
                &resultdir,
                sample_list.clone(),
                &format,
                bin_qualities.clone(),
                is_paired,
                assembler.clone(),
                id,
            )
            .map_err(|e| {
                eprintln!("Error processing bin {:?}: {}", component, e);
                e
            })
        })
        .expect("Error during processing components");
    });
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
    sample_list: Vec<String>,
    format: &String,
    bin_qualities: HashMap<String, utility::BinQuality>,
    is_paired: bool,
    assembler:String,
    id: usize,
) -> std::io::Result<()> {

    // eg: comp = {"binname_S1", "binname_S2"}
    debug!("{} id with length {} for component {:?}", id, component.len(), component);
    if component.len() == 1 {
        debug!("Single component. Process it further {:?}", component);
        let binname = component.into_iter().next().expect("The component is empty.");
        if let Some(quality) = bin_qualities.get(&binname) {
            if quality.completeness > 50.0 {
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
    
    if utility::check_high_quality_bin(&component, &bin_qualities, bindir, resultdir, id, &format) {
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

    utility::combine_fastabins(
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

    // eg. checkm2_subsetpath = <bindir>/0_checkm2_results/
    // let checkm2_subsetpath: PathBuf = 
    //     bindir.join(format!("{}_checkm2_results", id.to_string()));
    // let checkm2_subsetqualities = 
    //     utility::assess_bins(&selected_binset_path, 
    //     &checkm2_subsetpath, 8, "fasta").map_err(|e| {
    //         eprintln!(
    //             "Error in assessing combined bins for component {}: {}",
    //             id, e
    //         );
    //         e
    //     })?;

    // let mergebin_qualities = match utility::parse_bins_quality(
    //     &checkm2_subsetqualities
    // ) {
    //     Ok(quality) => quality,
    //     Err(_) => {
    //         eprintln!(
    //             "Failed to parse quality for mergedbin {:?}. Skipping bin.",
    //             &selected_binset_path
    //         );
            
    //         let _ = utility::select_highcompletebin(
    //             &component,
    //             &bin_qualities,
    //             &bindir,
    //             resultdir,
    //         );
    //         return Ok(());
    //     }
    // };

    // Check quality of combined bin and process only if it is high purity.
    // Otherwise choose the best bin
    // if let Some(bin_quality) = mergebin_qualities.get("combined") {
        
    //     // combined bin has low contamination, it will be processed further
    //     // TODO: If combined bin completeness is lower than unmerged bins, decision should be taken
    //     if bin_quality.contamination < 50f64 {
    //         // eg: from <bindir>/0_combined/combined.fasta
    //         // to <bindir>/0_combined/0_combined.fasta
    //         let _ = rename(selected_binset_path
    //             .join("combined.fasta"), 
    //             selected_binset_path.join(format!("{}_combined.fasta",id.to_string())));
    //     }
    // }
    
    // if !selected_binset_path.join(format!("{}_combined.fasta",id.to_string())).exists() {
    //     info!("For bin {}, no merged bin passes quality criteria to proceed to assembly", id);
    //     let _ = utility::select_highcompletebin(
    //         &component,
    //         &bin_qualities,
    //         &bindir,
    //         resultdir,
    //     );
    //     return Ok(());
    // }
    debug!("input to fetch combined contig ids {}", &selected_binset_path.join("combined.fasta").to_string_lossy());
    let mut all_enriched_scaffolds = HashSet::new();
    if !gfa_flag {
        all_enriched_scaffolds = gfaparser::read_fasta(
            &selected_binset_path.join("combined.fasta").to_string_lossy()
        )?;
    } else {
        let mut create_new = true;
        for samplebin in component.clone() {
            debug!("going to read gfa processing for {}", samplebin);
            let sample = sample_list
                .iter()
                .find(|sample| samplebin.contains(*sample))
                .cloned()
                .ok_or_else(|| {
                    error!("Bin in the component {} does not have a matching sample ID", samplebin);
                    std::io::Error::new(std::io::ErrorKind::NotFound, "Sample ID not found")
                })?;

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
        let sample = sample_list
            .iter()
            .find(|sample| samplebin.contains(*sample))
            .cloned()
            .ok_or_else(|| {
                error!("Bin in the component {} does not have a matching sample ID", samplebin);
                std::io::Error::new(std::io::ErrorKind::NotFound, "Sample ID not found")
            })?;
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

    let _ = utility::run_reassembly(
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
    );
    info!("Assembly step is gone through {}", id);
        // clean folders
        // let _ = fs::remove_dir_all(reassembly_outputdir);
        // let pattern = &format!("{}_combined*", i.to_string());
        // let _ = utility::remove_files_matching_pattern(
        //     &mergedbinpath,
        //     pattern);

    // clean folders
    // let pattern = &format!("*.{}", format);
    // let _ = utility::remove_files_matching_pattern(&binspecificdir, pattern);
    // let pattern = &format!("*.fastq");
    // let _ = utility::remove_files_matching_pattern(&binspecificdir, pattern);
    // let pattern = &format!("*_scaffolds*");
    // let _ = utility::remove_files_matching_pattern(&binspecificdir, pattern);
    // let _ = fs::remove_dir_all(checkm2_outputpath);
    // let _ = fs::remove_dir_all(binspecificdir);

    Ok(())
}
