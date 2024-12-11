use clap::{Arg, Command};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File, rename};
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use csv::ReaderBuilder;
use bio::io::fasta;
use petgraph::graph::Graph;
use petgraph::visit::Dfs;
use petgraph::Undirected;
use std::process::{Command as ProcessCommand, Stdio, exit};

mod parser;

use parser::parse_gfa_fastq;

#[derive(Debug)]
struct BinQuality {
    completeness: f64,
    contamination: f64,
}

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
    let bindir = validate_path(
        matches.
        get_one::<PathBuf>("bindir"), "bindir");
    let gfadir = validate_path(
        matches
        .get_one::<PathBuf>("gfadir"), "gfadir");
    let assemblydir = validate_path(
        matches
        .get_one::<PathBuf>("assemblydir"), "assemblydir");
    let mapdir = validate_path(
        matches
        .get_one::<PathBuf>("mapdir"), "mapdir");
    let readdir = validate_path(
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


    let binfiles = get_binfiles(bindir,format)?;

    for bin in binfiles {
        let binspecificdir = bindir.join(
            bin.file_name()
            .unwrap_or_default()
            .to_str()
            .unwrap_or_default())
            .with_extension("");
        let _ = splitbysampleid(&bin, &binspecificdir);
        let mergedbinpath: PathBuf = binspecificdir
            .join("mergedbins");
        // if mergedbinpath.exists() {
        //     fs::remove_dir_all(&mergedbinpath)?;
        // }
        // fs::create_dir(&mergedbinpath)?;
            
        let qualitybinpath: PathBuf = binspecificdir
            .join("qualitybins");
        if qualitybinpath.exists() {
            fs::remove_dir_all(&qualitybinpath)?;
        }
        fs::create_dir(&qualitybinpath)?;
        
        // *** Obtain quality of bins ***
        let checkm2_outputpath: PathBuf = binspecificdir
            .join("checkm2_results");
        
        let checkm2_qualities = assess_bins(
            &binspecificdir,
            &checkm2_outputpath,
            threads,
            format)?;
        let bin_qualities = parse_bins_quality(
            bindir,
            &qualitybinpath,
            &checkm2_qualities,
            format)?;
                
                
        // *** Find connected samples ***
        let bin_names: HashSet<String> = bin_qualities.keys().cloned().collect();
        combine_fastabins(&binspecificdir, &bin_names, format, &binspecificdir)?;
        let cluster_output = binspecificdir.join("linclust");
        
        let graph = find_overlappingbins(
            &&binspecificdir.join(format!("combined.{}", format)),
            cluster_output,
            min_overlaplen)?;
        let connected_samples = get_connected_samples(&graph);

        for (i, comp) in connected_samples.iter().enumerate() {

            // *** check quality of the components if merged ***
            let selected_binset_path = 
                mergedbinpath.join(format!("{}_combined", i.to_string()));
            if selected_binset_path.exists() {
                fs::remove_dir_all(&selected_binset_path)?;
            }
            fs::create_dir(&selected_binset_path)?;
            combine_fastabins(
                &binspecificdir,
                comp,
                format,
                &selected_binset_path)?;
            let checkm2_subsetpath: PathBuf = 
                mergedbinpath.join(format!("{}_checkm2_results", i.to_string()));
            let checkm2_qualities = 
                assess_bins(&selected_binset_path, 
                    &checkm2_subsetpath, threads, format)?;
            let mergebin_qualities = parse_bins_quality(
                &selected_binset_path,
                &checkm2_subsetpath,
                &checkm2_qualities,
                format
            )?;

            if let Some(bin_quality) = mergebin_qualities.get("combined") {
                // combined bins has low contamination, will be processed further
                if bin_quality.contamination < 10f64 {
                    let _ = rename(selected_binset_path.join(format!("combined.{}",format)), 
                    mergedbinpath.join(format!("{}_combined.{}",i.to_string(), format)));
                    let _ = fs::remove_dir_all(&checkm2_subsetpath);
                    let _ = fs::remove_dir_all(&selected_binset_path);
                } else {
                    let _ = select_highcompletebin(
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
            
            let mut create_new = true;
            // TODO: As of now, it collects information from only source samples of bins.
            // All samples would be better option and needs to be tested.
            for sample in comp {
                let gfa_path = gfadir.join(format!("{}.gfa", sample));
                let gfa_file = path_to_str(&gfa_path);

                let subbin_path = binspecificdir.join(format!("{}.{}", sample, format));
                let subbin_file = path_to_str(&subbin_path);
                
                let assembly_path = assemblydir.join(format!("{}.{}", sample, format));
                let assembly_file = path_to_str(&assembly_path);
                
                let mapid_path = mapdir.join(format!("{}_mapids", sample));
                let mapid_file = path_to_str(&mapid_path);

                let read_path = readdir.join(format!("{}.fastq", sample));
                let read_file = path_to_str(&read_path);
                
                let _ = parse_gfa_fastq(
                    gfa_file,
                    subbin_file,
                    assembly_file,
                    mapid_file,
                    read_file,
                    mergedbinpath
                    .join(
                    format!("{}_combined.{}",
                    i.to_string(), format)),
                    create_new
                );
                create_new = false;
            }
                
            let reassembly_outputdir = mergedbinpath.join(format!("{}_assembly",i.to_string()));
            let status = run_reassembly(
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
                    let _ = select_highcompletebin(
                        comp,
                        &bin_qualities,
                        &binspecificdir,
                        &mergedbinpath,
                        i,
                        format);
                        println!("Selected the best bin among samples based on completeness");
                }
                Err(_) => todo!(),
            }

            // clean folders
            
        }
        println!("{:?} successfully merged bins for input ", bindir);
    }

    Ok(())
}

fn validate_path<'a>(path: Option<&'a PathBuf>, name: &'a str) -> &'a PathBuf {
    let path = path.expect(&format!("{} path is required", name));
    if !path.exists() {
        eprintln!("Error: The specified path for {} does not exist", name);
        std::process::exit(1);
    }
    path
}

fn path_to_str(path: &PathBuf) -> &str {
    path.to_str().expect("Failed to convert PathBuf to &str")
}

fn splitbysampleid(
    bin: &PathBuf,
    binspecificdir: &PathBuf,
) -> io::Result<()>{

    if !binspecificdir.exists() {
        fs::create_dir_all(&binspecificdir)?;
    }
    // Open the input file
    let reader = BufReader::new(File::open(&bin)?);

    // Create a HashMap to store writers for each sample ID
    let mut writers: HashMap<String, File> = HashMap::new();

    let mut current_sample_id = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // Parse the sample ID from the header line
            if let Some(idx) = line.find('C') {
                current_sample_id = line[1..idx].to_string(); // Exclude the '>'

                // Ensure we have a writer for this sample ID
                if !writers.contains_key(&current_sample_id) {
                    let output_filename = binspecificdir.join(format!("{}.fasta", current_sample_id));
                    let output_file = File::create(output_filename)?;
                    writers.insert(current_sample_id.clone(), output_file);
                }
            } else {
                eprintln!("Warning: Could not find 'C' in header: {}", line);
            }
        }

        // Write the line to the appropriate file
        if let Some(writer) = writers.get_mut(&current_sample_id) {
            writeln!(writer, "{}", line)?;
        }
    }

    Ok(())
}

fn get_binfiles(dir: &Path, extension: &str) -> io::Result<Vec<PathBuf>> {
    let mut files = Vec::new();

    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();

        if path.is_file() {
            let filename = path.file_name().unwrap_or_default().to_str().unwrap_or_default();
            if filename.contains("_all_seqs") || filename.contains("rep_seq") || filename.contains("combined") {
                continue;
            }

            if let Some(ext) = path.extension() {
                if ext == extension {
                    files.push(path);
                }
            }
        }
    }

    Ok(files)
}

fn assess_bins(
    bindir: &PathBuf,
    bincheckm2: &PathBuf,
    threads: usize,
    format: &str,
) -> Result<PathBuf, io::Error> {
    let checkm2_qualities = Path::new(bincheckm2).join("quality_report.tsv");
        
    // Run CheckM2 run
    if !checkm2_qualities.exists() {
        println!("quality_report.tsv not found. Running checkm2...");
        
        let output = ProcessCommand::new("checkm2")
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
            .status()?;

        if !output.success() {
            eprintln!("Error: checkm2 command failed");
            exit(1);
        }
    }
    Ok(checkm2_qualities)
}

fn parse_bins_quality(
    inputdir: &PathBuf,
    movebinpath: &PathBuf,
    checkm2_qualities: &PathBuf,
    format: &str,
) -> io::Result<HashMap<String, BinQuality>> {


    // read checkm2 output file
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(File::open(checkm2_qualities)?);
    let mut bin_qualities: HashMap<String, BinQuality> = HashMap::new();

    for result in rdr.records() {
        let record = result?; // Get the record from the CSV
        if record.len() < 3 {
            eprintln!("Skipping invalid record: {:?}", record);
            continue; // Skip records that do not have enough columns
        }
        let bin_name = record[0].to_string();
        let completeness: f64 = record[1].parse().unwrap_or(0.0);
        let contamination: f64 = record[2].parse().unwrap_or(0.0);

        // Filter bins by contamination
        if contamination < 5.0f64 {
            if completeness > 95f64 {
                let bin_file_path = inputdir.join(format!("{}.{}", bin_name, format));
                let destination_path = movebinpath.join(format!("{}.{}", bin_name, format));
                fs::copy(&bin_file_path, &destination_path)?;
                println!("Copied existing high quality (96% comp, <5% cont) bin '{}' to '{}'", 
                    bin_name, movebinpath.display());
            } else {
                bin_qualities.insert(bin_name, 
                BinQuality{ completeness, contamination });
            }
        }
    }
    Ok(bin_qualities)
}

fn combine_fastabins(
    inputdir: &Path,
    bin_names: &HashSet<String>,
    format: &str,
    combined_bins: &Path,
) -> io::Result<()> {
    // Combine bins fasta into a single file
    let mut output_writer = File::create(combined_bins.join(format!("combined.{}", format)))?;
    for bin_name in bin_names {
        let bin_file_path = inputdir.join(format!("{}.{}", bin_name, format));

        if bin_file_path.exists() {
            let bin_file = File::open(&bin_file_path)?;
            let reader = fasta::Reader::new(bin_file);

            for record in reader.records() {
                let record = record?;  // Get the record
                writeln!(output_writer, ">{}", format!("{}",record.id()))?;
                writeln!(output_writer, "{}", String::from_utf8_lossy(record.seq()))?;
            }
        } else {
            eprintln!("Warning: File for bin '{}' does not exist at {:?}", bin_name, bin_file_path);
        }
    }
    Ok(())
}

fn find_overlappingbins(
    bins: &PathBuf,
    cluster_output: PathBuf,
    min_overlaplen: usize,
) -> io::Result<Graph<String, (), Undirected>> {
    let output = ProcessCommand::new("mmseqs")
        .arg("easy-linclust")
        .arg(bins)
        .arg(cluster_output.clone())
        .arg("tmp")
        .arg("--min-seq-id")
        .arg("1.0")
        .arg("--min-aln-len")
        .arg(min_overlaplen.to_string())
        .arg("--kmer-per-seq-scale")
        .arg("0.3")
        .stdout(Stdio::null())
        // .stderr(Stdio::null())
        .status()?; // Execute the command and get the status

    if !output.success() {
        eprintln!("Error: mmseqs easy-linclust failed");
        exit(1); // Exit if the command fails
    }

    let mut graph = Graph::<String, (), Undirected>::default();
    let mut bin_name_to_node = HashMap::new();
    let file = File::open(cluster_output.display().to_string() + "_cluster.tsv")?;
    let reader = io::BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let columns: Vec<&str> = line.split('\t').collect();

        // Assuming the bin names are in the first two columns (contig1, contig2)
        if columns.len() >= 2 {
            let bin1 = columns[0]
                .to_string()
                .split("C").next()
                .unwrap_or("").to_string();
            let bin2 = columns[1]
                .to_string()
                .split("C").next()
                .unwrap_or("").to_string();

            if bin1 != bin2 {
                let node1 = *bin_name_to_node
                    .entry(bin1.clone())
                    .or_insert_with(|| graph.add_node(bin1.clone()));
                let node2 = *bin_name_to_node
                    .entry(bin2.clone())
                    .or_insert_with(|| graph.add_node(bin2.clone()));

                graph.add_edge(node1, node2, ());
            }
        }
    }
    let files_to_remove = vec![
        cluster_output.display().to_string() + "_all_seqs.fasta",
        cluster_output.display().to_string() + "_rep_seq.fasta",
        cluster_output.display().to_string() + "_cluster.tsv"];
    for file_name in files_to_remove {
        let _ = fs::remove_file(file_name);
    }
   Ok(graph)
}

fn get_connected_samples(
    graph: &Graph<String, (), Undirected>
) -> Vec<HashSet<String>> {
    let mut visited = HashSet::new();
    let mut connected_samples = Vec::new();

    for node_index in graph.node_indices() {
        if !visited.contains(&node_index) {
            // Start a new component
            let mut component = HashSet::new();
            let mut dfs = Dfs::new(&graph, node_index);

            while let Some(nx) = dfs.next(&graph) {
                if visited.insert(nx) {
                    let node_name = graph[nx].clone();
                    component.insert(node_name);
                }
            }

            connected_samples.push(component);
        }
    }

    connected_samples
}

fn select_highcompletebin(
    comp: &HashSet<String>,
    bin_qualities: &HashMap<String, BinQuality>,
    binspecificdir: &PathBuf,
    mergedbinpath: &PathBuf,
    _i: usize,
    format: &str
) -> io::Result<()> {
    let highest_completebin = comp
        .iter()
        .filter_map(|bin
        | bin_qualities.get(bin)
        .map(|quality
        | (bin, quality.completeness)))
        .max_by(|(_, completeness1),(_, completeness2)
        | completeness1
        .partial_cmp(completeness2)
        .unwrap_or(std::cmp::Ordering::Equal))
        .map(|(bin, _)| bin.clone());
    
    if let Some(bin) = highest_completebin {
        let bin_path = binspecificdir.join(format!("{}.{}", bin, format));
        println!("{:?} {:?} highest bin complete", bin_path, format);
        rename(bin_path, mergedbinpath.join(format!("{}_final.{}",_i.to_string(), format)))?;
    } else {
        eprintln!("Error: No bin found with highest completeness.");
    }
    Ok(())
}

fn run_reassembly(
    readfile: &[PathBuf],
    binfile: &PathBuf,
    outputdir: &PathBuf,
    interleaved: bool,
    threads: usize,
) -> io::Result<bool> {
    
    if readfile.is_empty() {
        eprintln!("Error: No read files provided.");
        exit(1);
    }
    // Run SPAdes assembler
    let mut output = ProcessCommand::new("spades.py");
        output.arg("--trusted-contigs")
        .arg(binfile)
        .arg("--only-assembler")
        .arg("--careful")
        .arg("-o")
        .arg(outputdir)
        .arg("-t")
        .arg(threads.to_string())
        .arg("-m")
        .arg(128.to_string())
        .stdout(Stdio::null());
    if interleaved {
        output.arg("--12").arg(&readfile[0]);
    } else {
        output.arg("-1").arg(&readfile[0]).arg("-2").arg(&readfile[1]);
    } 

    match output.status() {
        Ok(status) if status.success() => {
            println!("SPAdes completed successfully.");
            Ok(true)
        }
        Ok(_) => {
            eprintln!("This is likely due to small input bin and insufficient k-mer counts from readset.");
            Ok(false)
        }
        Err(e) => {
            eprintln!("Error: Failed to execute SPAdes command - {}", e);
            Err(e)
        }
    }
}