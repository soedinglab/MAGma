use csv::ReaderBuilder;
use bio::io::fasta;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::fs::{self, File, rename};
use std::io::{self, BufRead, BufReader, Write};
use petgraph::graph::Graph;
use petgraph::visit::Dfs;
use petgraph::Undirected;
use std::process::{Command as ProcessCommand, Stdio, exit};
use glob::glob;

pub struct BinQuality {
    pub completeness: f64,
    pub contamination: f64,
}

pub fn validate_path<'a>(path: Option<&'a PathBuf>, name: &'a str, suffix: &str) -> &'a PathBuf {
    let path = path.expect(&format!("{} path is required", name));
    
    if !path.exists() {
        eprintln!("Error: The specified path for {} does not exist", name);
        std::process::exit(1);
    }

    if !path.is_dir() {
        eprintln!("Error: The specified path for {} is not a directory", name);
        exit(1);
    }

    let contains_files_with_suffix = fs::read_dir(path)
        .expect("Failed to read directory")
        .filter_map(|entry| entry.ok()) // Filter out invalid entries
        .any(|entry| {
            if let Some(file_name) = entry.file_name().to_str() {
                file_name.contains(suffix)
            } else {
                false
            }
        });

    if !contains_files_with_suffix {
        eprintln!(
            "Error: The directory for {} does not contain any files with the required extention/suffix '{}'",
            name, suffix
        );
        std::process::exit(1);
    }

    path
}

pub fn path_to_str(path: &PathBuf) -> &str {
    path.to_str().expect("Failed to convert PathBuf to &str")
}

pub fn check_paired_reads(directory: &PathBuf) -> bool {
    fs::read_dir(directory)
        .ok()
        .and_then(|entries
        | { entries
            .filter_map(|entry
            | entry.ok()?.file_name()
            .to_str().map(String::from))
            .find(|name
            | name.contains("_1") || name.contains("_2"))
    })
    .is_some()
}

pub fn find_file_with_extension(directory: &PathBuf, base_name: &str) -> PathBuf {
    let fastq = directory.join(format!("{}.fastq", base_name));
    if fastq.exists() {
        fastq
    } else {
        directory.join(format!("{}.fastq.gz", base_name))
    }
}

pub fn get_binfiles(dir: &Path, extension: &str) -> io::Result<Vec<PathBuf>> {
    let mut files = Vec::new();

    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();

        if path.is_file() {
            let filename = path.file_name()
                .unwrap_or_default()
                .to_str().unwrap_or_default();
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

pub fn splitbysampleid(
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
            current_sample_id = extract_sample_id(&line)?;
            ensure_writer(&current_sample_id, binspecificdir, &mut writers)?;
        }
        write_line_to_file(&current_sample_id, &line, &mut writers)?;
    }

    Ok(())
}

pub fn extract_sample_id(line: &str) -> io::Result<String> {
    if let Some(idx) = line.find('C') {
        Ok(line[1..idx].to_string()) // Exclude the '>'
    } else {
        eprintln!("Warning: Could not find 'C' in header: {}", line);
        Err(
        io::Error::new(
        io::ErrorKind::InvalidData
        ,"Invalid header format"))
    }
}


pub fn ensure_writer(
    sample_id: &str,
    binspecificdir: &Path,
    writers: &mut HashMap<String, File>,
) -> io::Result<()> {
    if !writers.contains_key(sample_id) {
        let output_filename = 
            binspecificdir
            .join(
            format!("{}.fasta"
            ,sample_id));
        let output_file =
            File::create(output_filename)?;
        writers.insert(
        sample_id
        .to_string()
        , output_file);
    }
    Ok(())
}

pub fn write_line_to_file(
    sample_id: &str,
    line: &str,
    writers: &mut HashMap<String, File>,
) -> io::Result<()> {
    if let Some(writer) =
        writers.get_mut(sample_id) {
            writeln!(writer, "{}", line)?;
    }
    Ok(())
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
                
                // eprintln!("NOTE: CheckM2 failed for {:?}! Check (-x) input format is incorrect."
                // ,bindir);
            }
            Err(e) => {
                eprintln!("Error: Failed to execute CheckM2 command - {}. Check if CheckM2 is executable currently", e);
            }
        }
    }    
    Ok(checkm2_qualities)
}

pub fn parse_bins_quality(
    inputdir: &PathBuf,
    movebinpath: &PathBuf,
    checkm2_qualities: &PathBuf,
    movebin_flag: bool,
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
        println!("bin_name {:?} completeness {:?} contamination {:?}", bin_name, completeness, contamination);
        // Filter bins by contamination
        if contamination < 5.0f64 {
            if completeness > 95f64 && movebin_flag {
                let bin_file_path = inputdir.join(format!("{}.fasta", bin_name));
                let destination_path = movebinpath.join(format!("{}.fasta", bin_name));
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

pub fn combine_fastabins(
    inputdir: &Path,
    bin_names: &HashSet<String>,
    combined_bins: &Path,
) -> io::Result<()> {
    // Combine bins fasta into a single file
    let mut output_writer = File::create(
        combined_bins
        .join(format!("combined.fasta")))?;
    for bin_name in bin_names {
        let bin_file_path = inputdir.join(format!("{}.fasta", bin_name));

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

pub fn find_overlappingbins(
    bins: &PathBuf,
    cluster_output: PathBuf,
    min_overlaplen: usize,
    threads: usize,
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
        .arg("--threads")
        .arg(threads.to_string())
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

pub fn get_connected_samples(
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

pub fn select_highcompletebin(
    comp: &HashSet<String>,
    bin_qualities: &HashMap<String, BinQuality>,
    binspecificdir: &PathBuf,
    outputpath: &PathBuf,
    _i: usize,
    bin_name: &str,
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
    
    if let Some(bin_sample) = highest_completebin {
        let bin_path = binspecificdir.join(format!("{}.fasta", bin_sample));
        println!("{:?} highest bin complete", bin_path);
        rename(bin_path, outputpath.join(format!("{}_{}_final.fasta",bin_name, _i.to_string())))?;
    } else {
        eprintln!("Error: No bin found with highest completeness.");
    }
    Ok(())
}

pub fn run_reassembly(
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
            eprintln!("NOTE: SPAdes failed for {:?}! This is likely due to small input bin and insufficient k-mer counts from readset."
                ,binfile.iter()
                .rev()
                .nth(2)  // nth(2) gives the third component (0-based index)
                .map(|bin| bin.to_string_lossy().to_string()).unwrap());
            Ok(false)
        }
        Err(e) => {
            eprintln!("Error: Failed to execute SPAdes command - {}", e);
            Err(e)
        }
    }
}

pub fn remove_files_matching_pattern(
    directory: &Path,
    pattern: &str
) -> io::Result<()> {
    let search_pattern = 
        directory
        .join(pattern)
        .to_string_lossy()
        .to_string();
    for entry in glob(&search_pattern)
        .map_err(|e
        | io::Error::new(
        io::ErrorKind::Other,
        e.to_string()))? {
        match entry {
            Ok(path) => {
                let _ = fs::remove_file(&path).is_err();
            }
            Err(e) => eprintln!("Failed to read matching file: {:?}", e),
        }
    }
    Ok(())
}