use clap::{Arg, Command};
use readfetch::readfetch;
use std::collections::{HashMap, HashSet};
use std::fs::{self, rename};
use std::io;
use std::path::PathBuf;
use rayon::prelude::*;
use log::{debug, info};

mod utility;
mod gfaparser;
mod readfetch;
use readfetch::fetch_fastqreads;
use gfaparser::parse_gfa_fastq;

fn main() -> io::Result<()> {

    env_logger::init();

    let matches = Command::new("mergebins")
        .version("1.0")
        .about("Merge redundant bins")
        .arg(
            Arg::new("bindir")
                .long("bindir")
                .short('b')
                .required(true)
                .help("Directory containing fasta files of bins")
                .value_parser(clap::value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("gfadir")
                .long("gfadir")
                .short('g')
                .required(true)
                .help("Directory containing gfa files for 
                    metagenomic samples (in gfa1.2 format)")
                .value_parser(clap::value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("assemblydir")
                .long("assemblydir")
                .short('a')
                .required(true)
                .help("Directory containing \
                    sample-wise assembly contigs file in fasta format")
                .value_parser(clap::value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("mapdir")
                .long("mapdir")
                .short('m')
                .required(true)
                .help("Directory containing mapids files \
                    derived from alignment sam/bam files \
                    (filename: <sample_id>_mapids)")
                .value_parser(clap::value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("readdir")
                .long("readdir")
                .short('r')
                .required(true)
                .help("Directory containing read files \
                    (filename: <sample_id>.fastq)")
                .value_parser(clap::value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("format")
                .long("format")
                .short('f')
                .default_value("fasta")
                .help("Bin file extension (fa/fas/fna/fasta), \
                    default=fasta"),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .short('t')
                .default_value("8")
                .help("Number of threads to use, default=8"),
        )
        .arg(
            Arg::new("min_overlaplen")
                .long("min_overlaplen")
                .short('l')
                .default_value("1000")
                .help("minimum overlap length, default=1000"),
        )
        .get_matches();

    
    // Parse arguments
    let format: &str = matches
    .get_one::<String>("format")
    .expect("Format is required");

    let bindir = utility::validate_path(
        matches.
        get_one::<PathBuf>("bindir"), "bindir", &format);
    let gfadir = utility::validate_path(
        matches
        .get_one::<PathBuf>("gfadir"), "gfadir", ".gfa");
    let assemblydir = utility::validate_path(
        matches
        .get_one::<PathBuf>("assemblydir"), "assemblydir", ".fasta");
    let mapdir = utility::validate_path(
        matches
        .get_one::<PathBuf>("mapdir"), "mapdir","_mapids");
    let readdir = utility::validate_path(
        matches
        .get_one::<PathBuf>("readdir"), "readdir", ".fastq");

    let threads: usize = matches
        .get_one::<String>("threads")
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(8);
    
    let min_overlaplen: usize = matches
        .get_one::<String>("min_overlaplen")
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(1000);


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
    
    let binfiles = utility::get_binfiles(bindir,format)?;

    if binfiles.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No bin files found. Please provide the correct format argument.",
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
    
    binfiles
        .par_iter()
        .filter(|bin| { // Optional pre-filtering step
            // Example filter: skip bins with invalid file names
            bin.file_name().is_some()
        })
        .try_for_each(|bin
        | process_bin(bin.clone(), 
            &bindir,
            &assemblydir,
            &gfadir,
            &mapdir,
            &readdir,
            &resultpath,
            &sample_list,
            threads,
            min_overlaplen,
            format,
            is_paired,
        ))?;

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
        &binspecificdir)?;

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
        );
        let _ = fs::remove_dir_all(binspecificdir);
        return Ok(());
    }
    debug!("At least two samples share overlapping region(s) {}", bin_name);
    let connected_samples = utility::get_connected_samples(&graph);
    
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
            &selected_binset_path)?;

        // eg. checkm2_subsetpath = <bindir>/bin_1/mergedbins/0_checkm2_results/
        let checkm2_subsetpath: PathBuf = 
            mergedbinpath.join(format!("{}_checkm2_results", i.to_string()));
        let checkm2_qualities = 
            utility::assess_bins(&selected_binset_path, 
            &checkm2_subsetpath, threads, format)?;
        let mergebin_qualities = utility::parse_bins_quality(
            &selected_binset_path,
            &checkm2_subsetpath,
            &checkm2_qualities,
            bin_name,
            false,
        )?;

        // Check quality of combined bin and process only if it is high purity.
        // Otherwise choose the best bin
        if let Some(bin_quality) = mergebin_qualities.get("combined") {
            
            // combined bin has low contamination, it will be processed further
            if bin_quality.contamination < 5f64 {
                // eg: from <bindir>/bin_1/mergedbins/0_combined/combined.fasta
                // to <bindir>/bin_1/mergedbins/0_combined.fasta
                let _ = rename(selected_binset_path
                    .join(format!("combined.{}",format)), 
                mergedbinpath.join(format!("{}_combined.{}",i.to_string(), format)));
            } else {
                // combined bin has high contamination. Choose the best
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
        // TODO: As of now, it collects information from only source samples of bins.
        // All samples would be better option and needs to be tested.
        // for sample in sample_list {
        for sample in comp {
            debug!("going to read gfa processing");
            // eg: <gfadir>/S1.gfa
            let gfa_path = gfadir.join(format!("{}.gfa", sample));
            let gfa_file = utility::path_to_str(&gfa_path);
            
            // eg: <bindir>/bin_1/S1.fasta
            let subbin_path = binspecificdir.join(format!("{}.fasta", sample));
            let subbin_file = utility::path_to_str(&subbin_path);
            
            // eg: <assemblydir>/S1.fasta. NOTE: contignames must contain sample id separated by 'C'
            let assembly_path = assemblydir.join(format!("{}.fasta", sample));
            let assembly_file = utility::path_to_str(&assembly_path);

            // eg: outputbin = <bindir>/bin_1/S1.fasta
            let enriched_scaffolds = parse_gfa_fastq(
                gfa_file,
                subbin_file,
                assembly_file,
                binspecificdir
                .join(
                format!("{}.fasta",
                sample)),
                create_new
            );
            for scaffold in enriched_scaffolds {
                // Check if the scaffold already contains the sample ID separated by 'C'.
                if !scaffold.contains(&format!("{}C", sample)) {
                    // Add sample ID to the scaffold if not present.
                    let modified_scaffold = format!("{}C{}", sample, scaffold);
                    all_enriched_scaffolds.insert(modified_scaffold);
                } else {
                    // If already present, keep the original scaffold.
                    all_enriched_scaffolds.insert(scaffold);
                }
            }
        }
        
        let create_new = true;
        for sample in sample_list {

            let mapid_path = mapdir.join(format!("{}_mapids", sample));
            let mapid_file = utility::path_to_str(&mapid_path);
            
            let (read_path1, read_path2): (PathBuf, Option<PathBuf>);
            let (read_file, read_file2): (&str, Option<&str>);
            
            if is_paired {
                read_path1 = utility::find_file_with_extension(&readdir, &format!("{}_1", sample));
                read_path2 = Some(utility::find_file_with_extension(&readdir, &format!("{}_2", sample)));
            
                read_file = utility::path_to_str(&read_path1);
                read_file2 = read_path2.as_ref().map(|path| utility::path_to_str(path));
            } else {
                read_path1 = utility::find_file_with_extension(&readdir, &sample);
                // read_path2 = None;
            
                read_file = utility::path_to_str(&read_path1);
                read_file2 = None;
            }

            let _ = fetch_fastqreads(
                &all_enriched_scaffolds,
                mapid_file,
                read_file,
                read_file2,
                binspecificdir
                .join(
                format!("{}.{}",
                sample, format)),
                is_paired,
                create_new
            );
            create_new = false;
        }
            
        let reassembly_outputdir = mergedbinpath.join(format!("{}_assembly",i.to_string()));
        let status = utility::run_reassembly(
            &[mergedbinpath.join(format!("{}_combined_enriched.fastq",i.to_string()))],
        &mergedbinpath.join(format!("{}_combined_enriched.{}",i.to_string(), format)),
        &reassembly_outputdir,
        true,
        threads,
        );

        match status {
            Ok(true) => {
                let _ = rename(reassembly_outputdir.join("scaffolds.fasta"), 
                resultpath.join(format!("{}_final.{}",i.to_string(), format)));
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
                let _ = rename(mergedbinpath.join(format!("{}_combined.{}",i.to_string(), format)), 
                resultpath.join(format!("{}_final.{}",i.to_string(), format)));
            }
            Err(_) => todo!(),
        }
        // clean folders
        let _ = fs::remove_dir_all(reassembly_outputdir);
        let pattern = &format!("{}_combined*", i.to_string());
        let _ = utility::remove_files_matching_pattern(
            &mergedbinpath,
            pattern);
    }
    // clean folders
    let pattern = &format!("*.{}", format);
    let _ = utility::remove_files_matching_pattern(&binspecificdir, pattern);
    let pattern = &format!("*.fastq");
    let _ = utility::remove_files_matching_pattern(&binspecificdir, pattern);
    let pattern = &format!("*_scaffolds*");
    let _ = utility::remove_files_matching_pattern(&binspecificdir, pattern);
    let _ = fs::remove_dir_all(checkm2_outputpath);
    let _ = fs::remove_dir_all(binspecificdir);

    Ok(())
}
