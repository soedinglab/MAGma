use std::collections::{HashMap, HashSet};
use std::fs::{self, rename};
use std::io::{self, stderr, Write};
use std::path::PathBuf;
use clap::Parser;
use rayon::prelude::*;
use rayon::{current_num_threads, ThreadPoolBuilder};
use log::{debug, error, info};
use num_cpus;

mod gfaparser;
mod readfetch;
mod utility;

// use index::indexfastqreads;
use readfetch::fetch_fastqreads;
use gfaparser::parse_gfa_fastq;

fn validate_paths(cli: &Cli) -> io::Result<(PathBuf, PathBuf, PathBuf, PathBuf, PathBuf)> {
    let bindir = utility::validate_path(Some(&cli.bindir), "bindir", &cli.format);
    let gfadir = utility::validate_path(Some(&cli.gfadir), "gfadir", ".gfa");
    let assemblydir = utility::validate_path(Some(&cli.assemblydir), "assemblydir", ".fasta");
    let mapdir = utility::validate_path(Some(&cli.mapdir), "mapdir", "_mapids");
    let readdir = utility::validate_path(Some(&cli.readdir), "readdir", ".fastq");

    Ok((bindir.to_path_buf(), gfadir.to_path_buf(), assemblydir.to_path_buf(), mapdir.to_path_buf(), readdir.to_path_buf()))
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
    gfadir: PathBuf,

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

    let sample_list = utility::get_gfa_file_names(&gfadir);
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
                let _ = utility::splitbysampleid(&bin, bin_name, &samplewisebinspath);
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
        "fasta")?;

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
        debug!("{:?} doesn't have high pure bins", bindir);
        return Ok(());
    }        

    let result = utility::calc_ani(&bindir, &bin_qualities, ani_cutoff);
    debug!("Graph was constructed for {:?}", bindir);
    let connected_bins: Vec<HashSet<String>> = match result {
        Ok(graph) => {
            debug!("Graph was constructed for {:?}", bindir);
            debug!("Graph {:?}", graph);
            // let bins = utility::get_connected_samples(&graph);
            let bins = utility::find_cliques(&graph);
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
                &assemblydir,
                &gfadir,
                &mapdir,
                &readdir,
                &resultdir,
                sample_list.clone(),
                bin_qualities.clone(),
                is_paired,
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
    assemblydir: &PathBuf,
    gfadir: &PathBuf,
    mapdir: &PathBuf,
    readdir: &PathBuf,
    resultdir: &PathBuf,
    sample_list: Vec<String>,
    bin_qualities: HashMap<String, utility::BinQuality>,
    is_paired: bool,
    id: usize,
) -> std::io::Result<()> {

    // eg: comp = {"binname_S1", "binname_S2"}
    debug!("{} id with length {} for component {:?}", id, component.len(), component);
    if component.len() == 1 {
        debug!("Single component. Process it further {:?}", component);
        let binname = component.into_iter().next().expect("The component is empty.");
        if let Some(quality) = bin_qualities.get(&binname) {
            if quality.completeness > 50.0 {
                let bin_path = bindir.join(format!("{}.fasta", binname));
                let final_path = resultdir.join(format!("{}.fasta", binname));
                debug!("Copying {:?} to {:?} for component {}", bin_path, final_path, id.to_string());
                
                if let Err(e) = fs::copy(&bin_path, &final_path) {
                    debug!("Failed to copy bin {}: {:?}", binname, e);
                }
            }
        }
        return Ok(());
    }
    
    

    if let Some((binname, best_quality)) = bin_qualities
        .iter()
        .filter(|(_, q)| q.completeness > 90.0) // Only consider completeness > 90
        .max_by(|a, b| a.1.completeness.partial_cmp(&b.1.completeness).unwrap()) 
    {
        info!(
            "Component {} has samplebin with highest completeness > 90: {} (Completeness: {})",
            id.to_string(), binname, best_quality.completeness
        );
        let bin_path = bindir.join(format!("{}.fasta", binname));
        let final_path = resultdir.join(format!("{}.fasta", binname));

        if let Err(e) = fs::copy(&bin_path, &final_path) {
            debug!("Failed to save bin {}: {:?}", binname, e);
        }

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
    &selected_binset_path).map_err(|e| {
        eprintln!(
            "Error in combining combined bins for component {}: {}",
            id, e
        );
        e
    })?;

    // eg. checkm2_subsetpath = <bindir>/0_checkm2_results/
    let checkm2_subsetpath: PathBuf = 
        bindir.join(format!("{}_checkm2_results", id.to_string()));
    let checkm2_subsetqualities = 
        utility::assess_bins(&selected_binset_path, 
        &checkm2_subsetpath, 8, "fasta").map_err(|e| {
            eprintln!(
                "Error in assessing combined bins for component {}: {}",
                id, e
            );
            e
        })?;

    let mergebin_qualities = match utility::parse_bins_quality(
        &checkm2_subsetqualities
    ) {
        Ok(quality) => quality,
        Err(_) => {
            eprintln!(
                "Failed to parse quality for mergedbin {:?}. Skipping bin.",
                &selected_binset_path
            );
            
            let _ = utility::select_highcompletebin(
                &component,
                &bin_qualities,
                &bindir,
                resultdir,
            );
            return Ok(());
        }
    };

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

    let _ = rename(selected_binset_path
        .join("combined.fasta"), 
        selected_binset_path.join(format!("{}_combined.fasta",id.to_string())));

    let mut create_new = true;
    let mut all_enriched_scaffolds = HashSet::new();
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
        let subbin_path = bindir.join(format!("{}.fasta", samplebin));
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
            .join(format!("{}_combined.fasta"
            ,id.to_string())),
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
    
    debug!("Fetched contigs from assembly file");
    create_new = true;
    debug!("all enriched scaffolds: {} {:?}", all_enriched_scaffolds.len(), all_enriched_scaffolds);
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
            .join(
            format!("{}_combined.fasta"
            ,id.to_string())),
            is_paired,
            create_new
        );
        create_new = false;
    }
    info!("Performing assembly for {}", id);
    let reassembly_outputdir = selected_binset_path.join(format!("{}_assembly",id.to_string()));
    let status = utility::run_reassembly(
        &[selected_binset_path.join(format!("{}_combined_enriched.fastq",id.to_string()))],
        &selected_binset_path.join(format!("{}_combined_enriched.fasta",id.to_string())),
        &reassembly_outputdir,
        true,
        8,
    );
    info!("Completed assembly for {}", id);
    debug!("status: {:?} {:?}", status, selected_binset_path.join("combined.fasta"));
    match status {
        Ok(true) => {
            let _ = utility::filterscaffold(&reassembly_outputdir.join("scaffolds.fasta"));
            let _ = rename(reassembly_outputdir.join("scaffolds_filtered.fasta"), 
            resultdir.join(format!("{}_merged_final.fasta", id.to_string())));
        }
        Ok(false) => {
            // select highest completeness bin
            let _ = utility::select_highcompletebin(
                &component,
                &bin_qualities,
                &bindir,
                &resultdir,
            );
            // fs::remove_dir_all(selected_binset_path)?;
            // fs::remove_dir_all(checkm2_subsetpath)?;
        }
        Err(_) => todo!(),
    }
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
