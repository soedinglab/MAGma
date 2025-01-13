use std::collections::HashSet;
use std::fs::{self, rename};
use std::io;
use std::path::PathBuf;
use rayon::prelude::*;
use rayon::{current_num_threads, ThreadPoolBuilder};
use log::{debug, info};

use crate::{index, utility, gfaparser, readfetch};
use index::indexfastqreads;
use readfetch::fetch_fastqreads;
use gfaparser::parse_gfa_fastq;


pub fn merge(
    bindir: PathBuf,
    gfadir: PathBuf,
    assemblydir: PathBuf,
    mapdir: PathBuf,
    readdir: PathBuf,
    format: String,
    threads: usize,
    min_overlaplen: usize,
) -> io::Result<()> {

    env_logger::init();
    
    // Parse arguments
    let bindir = utility::validate_path(Some(&bindir), "bindir", &format);
    let gfadir = utility::validate_path(Some(&gfadir), "gfadir", ".gfa");
    let assemblydir = utility::validate_path(Some(&assemblydir), "assemblydir", ".fasta");
    let mapdir = utility::validate_path(Some(&mapdir), "mapdir", "_mapids");
    let readdir = utility::validate_path(Some(&readdir), "readdir", ".fastq");

    println!("Starting merge process with the following parameters:");
    println!("  bindir: {:?}", bindir);
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
    
    let binfiles = utility::get_binfiles(bindir,&format)?;

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
    let resultpath: PathBuf = bindir
            .join("mergedbins");
    if resultpath.exists() {
        fs::remove_dir_all(&resultpath)?;
    }
    fs::create_dir(&resultpath)?;
    
    info!("{} cores are used", current_num_threads());

    let fastq_pool = ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .expect("Failed to build FASTQ thread pool");

    fastq_pool.install(|| {
        let _ = indexfastqreads(readdir, bindir); // Adjust arguments as needed
    });

    let pool = ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Failed to build thread pool");

    pool.install(|| {
        binfiles
        .par_iter()
        .filter(|bin| {
            // Optional pre-filtering step
            bin.file_name().is_some()
        })
        .try_for_each(|bin| {
            process_bin(
                bin.clone(),
                &bindir,
                &assemblydir,
                &gfadir,
                &mapdir,
                &readdir,
                &resultpath,
                &sample_list,
                threads,
                min_overlaplen,
                &format,
                is_paired,
            )
            .map_err(|e| {
                eprintln!(
                    "Error processing bin {:?}: {}",
                    bin.file_name().unwrap_or_default(),
                    e
                );
                e
            })
        })
        .expect("Error during processing");
    });
        
    Ok(())
    
}

fn process_bin(bin: PathBuf,
    bindir: &PathBuf,
    assemblydir: &PathBuf,
    gfadir: &PathBuf,
    mapdir: &PathBuf,
    readdir: &PathBuf,
    resultpath: &PathBuf,
    sample_list: &Vec<String>,
    threads: usize,
    min_overlaplen: usize,
    format: &str,
    is_paired: bool,
) -> std::io::Result<()> {

    // eg: <bindir>/bin_{1-N}.fasta
    // eg: bin_name = bin_1
    let bin_name = bin.file_name()
    .unwrap_or_default()
    .to_str()
    .unwrap_or_default().split('.').next().unwrap();
    
    // eg: binspecificdir = <bindir>/bin_1/
    let binspecificdir = bindir.join(
        bin_name)
        .with_extension("");
    let _ = utility::splitbysampleid(&bin, &binspecificdir);

    // Obtain quality of bins
    // eg: checkm2_outputpath = <bindir>/bin_1/checkm2_results/
    let checkm2_outputpath: PathBuf = binspecificdir
        .join("checkm2_results");
    
    let checkm2_qualities = utility::assess_bins(
        &binspecificdir,
        &checkm2_outputpath,
        threads,
        "fasta")?;

    let bin_qualities = utility::parse_bins_quality(
        &binspecificdir,
        &resultpath,
        &checkm2_qualities,
        bin_name,
        true)?;

    // None of the bins are pure
    if bin_qualities.len() == 0 {
        debug!("{:?} doesn't have high pure bins", bin_name);
        let _ = fs::remove_dir_all(binspecificdir);
        return Ok(());
    }        
        
    // Find connected samples
    let bin_samplenames: HashSet<String> = bin_qualities.keys().cloned().collect();
    debug!("{} bin_samplenames {:?}", bin_name, bin_samplenames);

    // if only single sample bin exist
    if bin_samplenames.len() == 1 {
        debug!("bin name for bin with a single sample {}", 
            bin_samplenames.iter().next().unwrap());
        debug!("result path {:?}",resultpath
        .join(format!("{}_final.fasta",bin_name)));

        // eg: from <bindir>/bin_1/S1.fasta to <bindir>/mergedbins/bin_1_final.fasta
        let _ = rename(binspecificdir
            .join(format!("{}.fasta",
            bin_samplenames.iter().next().unwrap())), 
            resultpath
            .join(format!("{}_final.fasta",bin_name)));
        let _ = fs::remove_dir_all(binspecificdir);
        return Ok(());
    }
    debug!("Multiple pure sample bins exist for bin {}", bin_name);
    utility::combine_fastabins(
        &binspecificdir, 
        &bin_samplenames,
        &binspecificdir).map_err(|e| {
            eprintln!(
                "Error in assess_bins for bin {}: {}",
                bin_name, e
            );
            e
        })?;

    // eg: cluster_output= <bindir>/bin_1/linclust
    let cluster_output = binspecificdir.join("linclust");
    debug!("cluster_output {:?} for {}", cluster_output, bin_name);
    
    // eg: combined.fasta = <bindir>/bin_1/combined.fasta
    let graph = utility::find_overlappingbins(
        &binspecificdir.join(format!("combined.fasta")),
        cluster_output,
        min_overlaplen,
        threads)?;
    debug!("Graph was constructed for {}", bin_name);
    // if no overlapping bins identified, 
    // only select the best bin by completeness.
    // Already bins were filtered by purity.
    if graph.node_count() == 0 {
        debug!("No overlapping bins have identified for {:?}. \
        Choosing the best bin by completeness", bin_name);
        let _ = utility::select_highcompletebin(
            &bin_samplenames,
            &bin_qualities,
            &binspecificdir,
            &resultpath,
            None,
            bin_name,
        ).map_err(|e| {
            eprintln!(
                "Error in assess_bins for bin {}: {}",
                bin_name, e
            );
            e
        })?;
        let _ = fs::remove_dir_all(binspecificdir);
        return Ok(());
    }
    debug!("At least two samples share overlapping region(s) {}", bin_name);
    let connected_samples = utility::get_connected_samples(&graph);
    
    if let Err(e) = utility::check_samplematch(&sample_list, &connected_samples) {
        eprintln!("Validation failed: {}", e);
    }

    // eg: <bindir>/bin_1/mergedbins/
    let mergedbinpath: PathBuf = binspecificdir
        .join("mergedbins");
    if mergedbinpath.exists() {
        fs::remove_dir_all(&mergedbinpath)?;
    }
    fs::create_dir(&mergedbinpath)?;

    debug!("Connected components {:?} for {}", connected_samples, bin_name);
    // eg: comp = {"S1", "S2"}
    for (i, comp) in connected_samples.iter().enumerate() {

        if comp.len() == 1 {
            debug!("Single component. Process it further {:?}", comp);
            // TODO:
            return Ok(());
        }

        // check quality of the components if merged
        // eg. selected_binset_path = <bindir>/bin_1/mergedbins/0_combined/
        let selected_binset_path = 
            mergedbinpath.join(format!("{}_combined", i.to_string()));
        if selected_binset_path.exists() {
            fs::remove_dir_all(&selected_binset_path)?;
        }
        fs::create_dir(&selected_binset_path)?;

        utility::combine_fastabins(
            &binspecificdir,
            comp,
            &selected_binset_path).map_err(|e| {
                eprintln!(
                    "Error in assess_bins for bin {}: {}",
                    bin_name, e
                );
                e
            })?;

        // eg. checkm2_subsetpath = <bindir>/bin_1/mergedbins/0_checkm2_results/
        let checkm2_subsetpath: PathBuf = 
            mergedbinpath.join(format!("{}_checkm2_results", i.to_string()));
        let checkm2_subsetqualities = 
            utility::assess_bins(&selected_binset_path, 
            &checkm2_subsetpath, threads, "fasta").map_err(|e| {
                eprintln!(
                    "Error in assess_bins for bin {}: {}",
                    bin_name, e
                );
                e
            })?;
        let mergebin_qualities = utility::parse_bins_quality(
            &selected_binset_path,
            &checkm2_subsetpath,
            &checkm2_subsetqualities,
            bin_name,
            false,
        ).map_err(|e| {
            eprintln!(
                "Error in assess_bins for bin {}: {}",
                bin_name, e
            );
            e
        })?;

        // Check quality of combined bin and process only if it is high purity.
        // Otherwise choose the best bin
        if let Some(bin_quality) = mergebin_qualities.get("combined") {
            
            // combined bin has low contamination, it will be processed further
            // TODO: If combined bin completeness is lower than unmerged bins, decision should be taken
            if bin_quality.contamination < 5f64 {
                // eg: from <bindir>/bin_1/mergedbins/0_combined/combined.fasta
                // to <bindir>/bin_1/mergedbins/0_combined.fasta
                let _ = rename(selected_binset_path
                    .join(format!("combined.{}",format)), 
                mergedbinpath.join(format!("{}_combined.fasta",i.to_string())));
            } else {
                // combined bin has high contamination. Choose the best
                // TODO: if unmerged bin is higher complete than merged one, choose that.
                let _ = utility::select_highcompletebin(
                    comp,
                    &bin_qualities,
                    &binspecificdir,
                    &mergedbinpath,
                    Some(i.to_string()),
                    bin_name,
                );
                return Ok(());
            }
            
        } else {
            debug!("Error in parsing merged bin checkm2 statistics");
        }
        
        let mut create_new = true;
        let mut all_enriched_scaffolds = HashSet::new();
        let db_path = bindir.join("reads.db");
        for sample in comp {
            debug!("going to read gfa processing");
            // eg: <gfadir>/S1.gfa
            let gfa_path = gfadir.join(format!("{}.gfa", sample));
            let gfa_file = utility::path_to_str(&gfa_path);
            
            // eg: <bindir>/bin_1/S1.fasta
            let subbin_path = binspecificdir.join(format!("{}.fasta", sample));
            let subbin_file = utility::path_to_str(&subbin_path);
            
            // eg: <assemblydir>/S1.fasta
            let assembly_path = assemblydir.join(format!("{}.fasta", sample));
            let assembly_file = utility::path_to_str(&assembly_path);

            // eg: outputbin = <bindir>/bin_1/S1.fasta
            let enriched_scaffolds = parse_gfa_fastq(
                gfa_file,
                subbin_file,
                assembly_file,
                mergedbinpath
                .join(format!("{}_combined.fasta"
                ,i.to_string())),
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
        // for sample in sample_list {
        // for sample in comp {

        //     debug!("processing sample file {} to get reads", sample);
        //     let mapid_path = mapdir.join(format!("{}_mapids", sample));
        //     let mapid_file = utility::path_to_str(&mapid_path);
            
        //     let (read_path1, read_path2): (PathBuf, Option<PathBuf>);
        //     let (read_file, read_file2): (&str, Option<&str>);
            
        //     if is_paired {
        //         read_path1 = utility::find_file_with_extension(&readdir, &format!("{}_1", sample));
        //         read_path2 = Some(utility::find_file_with_extension(&readdir, &format!("{}_2", sample)));
            
        //         read_file = utility::path_to_str(&read_path1);
        //         read_file2 = read_path2.as_ref().map(|path| utility::path_to_str(path));
        //     } else {
        //         read_path1 = utility::find_file_with_extension(&readdir, &sample);
            
        //         read_file = utility::path_to_str(&read_path1);
        //         read_file2 = None;
        //     }
            
        //     debug!("map id {} read file1 {:?} read file2 {:?}", mapid_file, read_file, read_file2);
        //     let _ = fetch_fastqreads(
        //         &all_enriched_scaffolds,
        //         mapid_file,
        //         read_file,
        //         read_file2,
        //         mergedbinpath
        //         .join(
        //         format!("{}_combined.fasta"
        //         ,i.to_string())),
        //         is_paired,
        //         create_new
        //     );
        //     create_new = false;
        // }

        for sample in comp {

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
                &db_path,
                mergedbinpath
                .join(
                format!("{}_combined.fasta"
                ,i.to_string())),
                create_new
            );
            create_new = false;
        }
        info!("Performing assembly for {}", bin_name);
        let reassembly_outputdir = mergedbinpath.join(format!("{}_assembly",i.to_string()));
        let status = utility::run_reassembly(
            &[mergedbinpath.join(format!("{}_combined_enriched.fastq",i.to_string()))],
            &mergedbinpath.join(format!("{}_combined_enriched.fasta",i.to_string())),
            &reassembly_outputdir,
            true,
            threads,
        );
        info!("Completed assembly for {}", bin_name);
        match status {
            Ok(true) => {
                let _ = utility::filterscaffold(&reassembly_outputdir.join("scaffolds.fasta"));
                let _ = rename(reassembly_outputdir.join("scaffolds_filtered.fasta"), 
                resultpath.join(format!("{}_final.fasta",i.to_string())));
            }
            Ok(false) => {
                // select highest completeness bin
                // let _ = select_highcompletebin(
                //     comp,
                //     &bin_qualities,
                //     &binspecificdir,
                //     &mergedbinpath,
                //     Some(i.to_string()),
                //     bin_name);
                //     println!("Selected the best bin among samples based on completeness");
                
                // output the combined bin
                let _ = rename(mergedbinpath.join(format!("{}_combined.fasta",i.to_string())), 
                resultpath.join(format!("{}_final.fasta",i.to_string())));
            }
            Err(_) => todo!(),
        }
        // clean folders
        // let _ = fs::remove_dir_all(reassembly_outputdir);
        // let pattern = &format!("{}_combined*", i.to_string());
        // let _ = utility::remove_files_matching_pattern(
        //     &mergedbinpath,
        //     pattern);
    }
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
