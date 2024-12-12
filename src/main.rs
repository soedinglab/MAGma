use clap::{Arg, Command};
use std::collections::HashSet;
use std::fs::{self, rename};
use std::io;
use std::path::PathBuf;

mod gfaparser;
use gfaparser::parse_gfa_fastq;
mod utility;

fn main() -> io::Result<()> {
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
                .help("Directory containing gfa files for metagenomic samples (in gfa1.2 format)")
                .value_parser(clap::value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("assemblydir")
                .long("assemblydir")
                .short('a')
                .required(true)
                .help("Directory containing sample-wise assembly contigs file in fasta format")
                .value_parser(clap::value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("mapdir")
                .long("mapdir")
                .short('m')
                .required(true)
                .help("Directory containing mapids files derived from alignment sam/bam files (filename: <sample_id>_mapids)")
                .value_parser(clap::value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("readdir")
                .long("readdir")
                .short('r')
                .required(true)
                .help("Directory containing read files (filename: <sample_id>.fastq)")
                .value_parser(clap::value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("format")
                .long("format")
                .short('f')
                .default_value("fasta")
                .help("Bin file extension (fa/fas/fna/fasta), default=fasta"),
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

    
    // *** Parse arguments *** 
    let bindir = utility::validate_path(
        matches.
        get_one::<PathBuf>("bindir"), "bindir");
    let gfadir = utility::validate_path(
        matches
        .get_one::<PathBuf>("gfadir"), "gfadir");
    let assemblydir = utility::validate_path(
        matches
        .get_one::<PathBuf>("assemblydir"), "assemblydir");
    let mapdir = utility::validate_path(
        matches
        .get_one::<PathBuf>("mapdir"), "mapdir");
    let readdir = utility::validate_path(
        matches
        .get_one::<PathBuf>("readdir"), "readdir");

    let format: &str = matches
        .get_one::<String>("format")
        .expect("Format is required");

    let threads: usize = matches
        .get_one::<String>("threads")
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(8);
    
    let min_overlaplen: usize = matches
        .get_one::<String>("min_overlaplen")
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(1000);


    let binfiles = utility::get_binfiles(bindir,format)?;

    for bin in binfiles {
        let binspecificdir = bindir.join(
            bin.file_name()
            .unwrap_or_default()
            .to_str()
            .unwrap_or_default())
            .with_extension("");
        let _ = utility::splitbysampleid(&bin, &binspecificdir);
        let mergedbinpath: PathBuf = binspecificdir
            .join("mergedbins");
        if mergedbinpath.exists() {
            fs::remove_dir_all(&mergedbinpath)?;
        }
        fs::create_dir(&mergedbinpath)?;
            
        let qualitybinpath: PathBuf = binspecificdir
            .join("qualitybins");
        if qualitybinpath.exists() {
            fs::remove_dir_all(&qualitybinpath)?;
        }
        fs::create_dir(&qualitybinpath)?;
        
        // *** Obtain quality of bins ***
        let checkm2_outputpath: PathBuf = binspecificdir
            .join("checkm2_results");
        
        let checkm2_qualities = utility::assess_bins(
            &binspecificdir,
            &checkm2_outputpath,
            threads,
            format)?;
        let bin_qualities = utility::parse_bins_quality(
            bindir,
            &qualitybinpath,
            &checkm2_qualities,
            format)?;
                
                
        // *** Find connected samples ***
        let bin_names: HashSet<String> = bin_qualities.keys().cloned().collect();
        utility::combine_fastabins(&binspecificdir, &bin_names, format, &binspecificdir)?;
        let cluster_output = binspecificdir.join("linclust");
        
        let graph = utility::find_overlappingbins(
            &&binspecificdir.join(format!("onlypurebins.{}", format)),
            cluster_output,
            min_overlaplen)?;
        let connected_samples = utility::get_connected_samples(&graph);
        println!("{:?}", connected_samples);

        for (i, comp) in connected_samples.iter().enumerate() {

            // *** check quality of the components if merged ***
            let selected_binset_path = 
                mergedbinpath.join(format!("{}_combined", i.to_string()));
            if selected_binset_path.exists() {
                fs::remove_dir_all(&selected_binset_path)?;
            }
            fs::create_dir(&selected_binset_path)?;
            utility::combine_fastabins(
                &binspecificdir,
                comp,
                format,
                &selected_binset_path)?;
            let checkm2_subsetpath: PathBuf = 
                mergedbinpath.join(format!("{}_checkm2_results", i.to_string()));
            let checkm2_qualities = 
                utility::assess_bins(&selected_binset_path, 
                    &checkm2_subsetpath, threads, format)?;
            let mergebin_qualities = utility::parse_bins_quality(
                &selected_binset_path,
                &checkm2_subsetpath,
                &checkm2_qualities,
                format
            )?;

            if let Some(bin_quality) = mergebin_qualities.get("combined") {
                // combined bin has low contamination, it will be processed further
                if bin_quality.contamination < 10f64 {
                    let _ = rename(selected_binset_path.join(format!("combined.{}",format)), 
                    mergedbinpath.join(format!("{}_combined.{}",i.to_string(), format)));
                } else {
                    let _ = utility::select_highcompletebin(
                        comp,
                        &bin_qualities,
                        &binspecificdir,
                        &mergedbinpath,
                        i,
                        format);
                        continue;
                }
                
            } else {
                println!("Error in parsing merged bin checkm2 statistics");
            }
            let _ = fs::remove_dir_all(&checkm2_subsetpath);
            let _ = fs::remove_dir_all(&selected_binset_path);
            
            let mut create_new = true;
            // TODO: As of now, it collects information from only source samples of bins.
            // All samples would be better option and needs to be tested.
            for sample in comp {
                let gfa_path = gfadir.join(format!("{}.gfa", sample));
                let gfa_file = utility::path_to_str(&gfa_path);

                let subbin_path = binspecificdir.join(format!("{}.{}", sample, format));
                let subbin_file = utility::path_to_str(&subbin_path);
                
                let assembly_path = assemblydir.join(format!("{}.{}", sample, format));
                let assembly_file = utility::path_to_str(&assembly_path);
                
                let mapid_path = mapdir.join(format!("{}_mapids", sample));
                let mapid_file = utility::path_to_str(&mapid_path);

                let read_path = readdir.join(format!("{}.fastq", sample));
                let read_file = utility::path_to_str(&read_path);
                
                let _ = parse_gfa_fastq(
                    gfa_file,
                    subbin_file,
                    assembly_file,
                    mapid_file,
                    read_file,
                    binspecificdir
                    .join(
                    format!("{}.{}",
                    sample, format)),
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
                    mergedbinpath.join(format!("{}_final.{}",i.to_string(), format)));
                }
                Ok(false) => {
                    // *** select highest completeness bin ***
                    // let _ = select_highcompletebin(
                    //     comp,
                    //     &bin_qualities,
                    //     &binspecificdir,
                    //     &mergedbinpath,
                    //     i,
                    //     format);
                    //     println!("Selected the best bin among samples based on completeness");
                    
                    // *** output the combined bin ***
                    let _ = rename(mergedbinpath.join(format!("{}_combined.{}",i.to_string(), format)), 
                    mergedbinpath.join(format!("{}_final.{}",i.to_string(), format)));
                }
                Err(_) => todo!(),
            }
            // clean folders
            let _ = fs::remove_dir_all(reassembly_outputdir);
            let _ = fs::remove_dir_all(mergedbinpath
                .join(
                format!("{}_combined*"
                ,i.to_string())));
            // let pattern = &format!("{}_combined*", i.to_string());
            // let _ = remove_files_matching_pattern(
            //     &mergedbinpath,
            //     pattern);
        }
        // clean folders
        let _ = fs::remove_file(binspecificdir.join(format!("*.{}",format)));
        let _ = fs::remove_file(binspecificdir.join("*.fastq"));
        let _ = fs::remove_file(binspecificdir.join("*_scaffolds*"));
        let _ = fs::remove_dir_all(checkm2_outputpath);
    }
    
    println!("{:?} successfully merged bins for input ", bindir);
    Ok(())
}
