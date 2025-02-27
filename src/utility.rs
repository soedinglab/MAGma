use csv::ReaderBuilder;
use bio::io::fasta;
use std::collections::{HashMap, HashSet, VecDeque};
use std::path::{Path, PathBuf};
use std::fs::{self, File, rename};
use std::io::{self, BufRead, BufReader, Write};
use petgraph::graph::Graph;
use petgraph::visit::Dfs;
use petgraph::{Undirected, prelude};
use std::process::{exit, Command as ProcessCommand, Stdio};
use log::{debug, info, error};
use glob::glob;

#[derive(Clone)]
pub struct BinQuality {
    pub completeness: f64,
    pub contamination: f64,
}

pub fn validate_path<'a>(path: Option<&'a PathBuf>, name: &'a str, suffix: &str) -> &'a PathBuf {
    let path = path.expect(&format!("{} path is required", name));
    
    if !path.exists() {
        error!("Error: The specified path for {} does not exist", name);
        std::process::exit(1);
    }

    if !path.is_dir() {
        error!("Error: The specified path for {} is not a directory", name);
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
        error!(
            "Error: The directory for {} does not contain any files with the required extention/suffix '{}'",
            name, suffix
        );
        std::process::exit(1);
    }

    path
}

pub fn path_to_str(path: &PathBuf) -> &str {
    path.to_str()
    .expect("Failed to convert PathBuf to &str")
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

pub fn get_gfa_file_names(dir: &Path) -> Vec<String> {
    // Try to read the directory entries
    match fs::read_dir(dir) {
        Ok(entries) => {
            // Filter and collect file names with .gfa extension
            entries
            .filter_map(|entry| entry.ok()) // Ignore errors when reading entries
            .filter(|entry| {
            entry
            .path()
            .extension()
            .and_then(|ext| ext.to_str()) // Get the extension as a string
            == Some("gfa")               // Compare the extension with "gfa"
            })
            .filter_map(|entry| {
                // Extract just the file name without the extension
            entry.path().file_stem().and_then(|name| name.to_str().map(|s| s.to_string()))
            })
            .collect() // Collect the file names into a Vec<String>
        }
        Err(_) => Vec::new(), // Return empty vector if directory reading fails
    }
}


pub fn splitbysampleid(
    bin: &PathBuf,
    bin_name: &str,
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
            ensure_writer(&current_sample_id, bin_name, binspecificdir, &mut writers)?;
        }
        write_line_to_file(&current_sample_id, &line, &mut writers)?;
    }
    debug!("Finished writing sample-wise bins for {:?}", bin);
    Ok(())
}

pub fn extract_sample_id(line: &str) -> io::Result<String> {
    if let Some(idx) = line.find('C') {
        Ok(line[1..idx].to_string()) // Exclude the '>'
    } else {
        error!("Warning: Could not find 'C' in header: {}", line);
        Err(
        io::Error::new(
        io::ErrorKind::InvalidData
        ,"Invalid header format"))
    }
}

pub fn ensure_writer(
    sample_id: &str,
    bin_name: &str,
    binspecificdir: &Path,
    writers: &mut HashMap<String, File>,
) -> io::Result<()> {
    if !writers.contains_key(sample_id) {
        let output_filename = 
            binspecificdir
            .join(
            format!("{}_{}.fasta"
            ,bin_name,sample_id));
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
                
                // error!("NOTE: CheckM2 failed for {:?}! Check (-x) input format is incorrect."
                // ,bindir);
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

pub fn combine_fastabins(
    inputdir: &Path,
    bin_samplenames: &HashSet<String>,
    combined_bins: &Path,
) -> io::Result<()> {
    // Combine bins fasta into a single file
    let mut output_writer = File::create(
        combined_bins
        .join(format!("combined.fasta")))?;
    for bin_samplename in bin_samplenames {
        let bin_file_path = inputdir.join(format!("{}.fasta", bin_samplename));
        debug!("Combine_fastabins: bin_file_path {:?} with {:?}", bin_file_path, combined_bins);
        if bin_file_path.exists() {
            let bin_file = File::open(&bin_file_path)?;
            let reader = fasta::Reader::new(bin_file);

            for record in reader.records() {
                let record = record?;  // Get the record
                writeln!(output_writer, ">{}", format!("{}",record.id()))?;
                writeln!(output_writer, "{}", String::from_utf8_lossy(record.seq()))?;
            }
        } else {
            error!("Warning: File for bin '{}' does not exist at {:?}", bin_samplename, bin_file_path);
        }
    }
    Ok(())
}

pub fn calc_ani(
    bins: &PathBuf,
    bin_qualities: &HashMap<String, BinQuality>,
    ani_cutoff: f32,
) -> io::Result<Graph<String, (), Undirected>> {
    let ani_output: PathBuf = bins.join("ani_edges");
    let bin_files: Vec<String> = glob(&format!("{}/*.fasta", bins.display()))
    .expect("Failed to read glob pattern")
    .filter_map(Result::ok)
    .map(|path| path.to_string_lossy().into_owned())
    .collect();

    if bin_files.is_empty() {
        error!("No fasta files found in {:?}", bins);
        return Err(io::Error::new(io::ErrorKind::NotFound, "No fasta files found"));
    }

    debug!("ani input to skani: {:?}", bin_files);
    let mut command = ProcessCommand::new("skani");
    command.arg("triangle");
    command.args(&bin_files);
    command.arg("-E");

    let output_file = File::create(&ani_output)?;
    command.stdout(Stdio::from(output_file));
    command.stderr(Stdio::inherit()); // Keep stderr visible

    let status = command.status()?;
    if !status.success() {
        return Err(io::Error::new(io::ErrorKind::Other, "skani triangle failed"));
    }
    
    let mut graph: Graph<String, (), Undirected> = Graph::default();
    let mut bin_name_to_node: HashMap<String, prelude::NodeIndex> = HashMap::new();
    let file: File = File::open(ani_output)?;
    let reader: BufReader<File> = io::BufReader::new(file);

    for line in reader.lines().skip(1) {
        let line: String = line?;
        let columns: Vec<&str> = line.split('\t').collect();

        let bin1 = Path::new(columns[0])
            .file_stem()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| columns[0].to_string());

        let bin2 = Path::new(columns[1])
            .file_stem()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| columns[1].to_string());
        debug!("bin1: {} bin2: {}", bin1, bin2);
        let ani: f32 = columns[2].parse().expect("Failed to parse ANI value as float from column 3");
        
        if !bin_qualities
            .get(&bin1)
            .map_or(false, |q| q.contamination < 5.0)
            || !bin_qualities
                .get(&bin2)
                .map_or(false, |q| q.contamination < 5.0)
        {
            continue;
        }
    
        let is_bin1_valid = bin_qualities
            .get(&bin1)
            .map_or(false, |q| q.completeness > 20.0);

        debug!("For bin1 {:?}, satisfy completeness {}", bin1, is_bin1_valid);
    
        let is_bin2_valid = bin_qualities
            .get(&bin2)
            .map_or(false, |q| q.completeness > 20.0);
        debug!("For bin2 {:?}, satisfy completeness {}", bin2, is_bin2_valid);
        
        let node1 = if is_bin1_valid {
            Some(*bin_name_to_node
                .entry(bin1.clone())
                .or_insert_with(|| graph.add_node(bin1.clone())))
        } else {
            None
        };
        
        let node2 = if is_bin2_valid {
            Some(*bin_name_to_node
                .entry(bin2.clone())
                .or_insert_with(|| graph.add_node(bin2.clone())))
        } else {
            None
        };
        
        // Add edge only if both nodes are valid, bins are different, and ANI > 99.0
        if let (Some(node1), Some(node2)) = (node1, node2) {
            if bin1 != bin2 && ani > ani_cutoff {
                graph.add_edge(node1, node2, ());
            }
        }
    }
   Ok(graph)
}

/// Find connected components that are cliques (fully connected subgraphs).
pub fn find_cliques(graph: &Graph<String, (), Undirected>) -> Vec<HashSet<String>> {
    let mut visited = HashSet::new();
    let mut clique_components = Vec::new();

    for node_index in graph.node_indices() {
        if !visited.contains(&node_index) {
            let mut component = HashSet::new();
            let mut queue = VecDeque::new();
            queue.push_back(node_index);

            while let Some(nx) = queue.pop_front() {
                if visited.insert(nx) {
                    let node_name = graph[nx].clone();
                    component.insert(node_name.clone());

                    // Add neighbors to queue
                    for neighbor in graph.neighbors(nx) {
                        if !visited.contains(&neighbor) {
                            queue.push_back(neighbor);
                        }
                    }
                }
            }

            // Check if the component is a clique (fully connected subgraph)
            if is_clique(graph, &component) {
                clique_components.push(component);
            }
        }
    }

    clique_components
}

/// Check if a given set of nodes form a clique (fully connected subgraph).
fn is_clique(graph: &Graph<String, (), Undirected>, component: &HashSet<String>) -> bool {
    let nodes: Vec<_> = component.iter().collect();
    
    // Every pair of nodes in the component should be directly connected
    for (i, node1) in nodes.iter().enumerate() {
        for node2 in nodes.iter().skip(i + 1) {
            let index1 = graph.node_indices().find(|&idx| &graph[idx] == *node1).unwrap();
            let index2 = graph.node_indices().find(|&idx| &graph[idx] == *node2).unwrap();

            if !graph.contains_edge(index1, index2) {
                return false; // Not a clique
            }
        }
    }
    true
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
    bin_samplenames: &HashSet<String>,
    bin_qualities: &HashMap<String, BinQuality>,
    bindir: &PathBuf,
    outputpath: &PathBuf,
) -> io::Result<()> {
    let highest_completebin = bin_samplenames
        .iter()
        .filter_map(|bin| bin_qualities.get(bin)
        .map(|quality| (bin, quality.completeness)))
        .max_by(|(_, completeness1), (_, completeness2)| {
            completeness1.partial_cmp(completeness2)
            .unwrap_or(std::cmp::Ordering::Equal)
        })
        .filter(|(_, max_completeness)| *max_completeness >= 50.0)
        .map(|(bin, _)| bin.clone());

    if let Some(sample_id) = highest_completebin {
        let bin_path = bindir.join(format!("{}.fasta", sample_id));
        debug!("Output path {:?}, Sample ID {:?}", outputpath, sample_id);
        
        fs::copy(bin_path, outputpath.join(format!("{}.fasta", sample_id)))?;
    }
    Ok(())

}


pub fn process_high_quality_bins(
    comp: &HashSet<String>,
    bin_qualities: &HashMap<String, BinQuality>,
    binspecificdir: &PathBuf,
    resultpath: &PathBuf,
    bin_name: &str,
) -> io::Result<()> {
    let comp_binqualities: HashMap<String, BinQuality> = bin_qualities
        .iter()
        .filter(|(bin, _)| comp.contains(*bin))
        .map(|(bin, quality)| (bin.to_string(), quality.clone()))
        .collect();

    // select_highcompletebin(
    //     comp,
    //     &comp_binqualities,
    //     binspecificdir,
    //     resultpath,
    //     sid,
    //     bin_name,
    // )
    // .map_err(|e| {
    //     eprintln!("Error in assess_bins for bin {}: {}", bin_name, e);
    //     e
    // })?;

    let _  = save_selectedbins(
        &comp_binqualities,
        binspecificdir,
        resultpath,
        bin_name
    );
    Ok(())
}

pub fn save_selectedbins(
    bin_qualities: &HashMap<String, BinQuality>,
    bindir: &PathBuf,
    outputpath: &PathBuf,
    bin_name: &str,
) -> io::Result<()> {
   
    let mut selected_bins: Vec<&String> = Vec::new();

    let thresholds = [
        (90.0, 5.0),  // First attempt: completeness > 90%, contamination < 5%
        (70.0, 10.0), // Second attempt: completeness > 70%, contamination < 10%
        // (50.0, 10.0), // Third attempt: completeness > 50%, contamination < 10%
    ];

    for &(completeness_thresh, contamination_thresh) in &thresholds {
        selected_bins = bin_qualities
            .iter()
            .filter(|(_, quality)| quality.completeness >= completeness_thresh 
                && quality.contamination < contamination_thresh)
            .map(|(bin, _)| bin)
            .collect();
        debug!("selected bins {:?} for {}", selected_bins, bin_name);
        if !selected_bins.is_empty() {
            break;
        }

        debug!(
            "No bins met the criteria (completeness > {}, contamination < {}) for {}",
            completeness_thresh, contamination_thresh, bin_name
        );
    }

    if selected_bins.is_empty() {
        debug!("No bins met the quality criteria for {}", bin_name);
        return Ok(()); // Exit early if no bins match
    }

    for bin in selected_bins {
        let bin_path = bindir.join(format!("{}.fasta", bin));
        let final_filename = format!("{}_{}.fasta", bin_name, bin);
        let final_path = outputpath.join(final_filename);

        if let Err(e) = rename(&bin_path, &final_path) {
            debug!("Failed to save bin {}: {:?}", bin, e);
        }
    }
    Ok(())
}


pub fn run_reassembly(
    readfile: &[PathBuf],
    binfile: &PathBuf,
    outputdir: &PathBuf,
    is_paired: bool,
    threads: usize,
) -> io::Result<bool> {
    
    if readfile.is_empty() {
        error!("Error: No read files provided.");
        exit(1);
    }
    // Run SPAdes assembler
    let mut output = ProcessCommand::new("spades.py");
        output.arg("--trusted-contigs")
        .arg(binfile)
        .arg("--only-assembler")
        // .arg("--careful")
        .arg("-o")
        .arg(outputdir)
        .arg("-t")
        .arg(threads.to_string())
        .arg("-m")
        .arg(128.to_string())
        .stdout(Stdio::null());
    if is_paired {
        output.arg("--12").arg(&readfile[0]);
    } else {
        output.arg("-1").arg(&readfile[0]).arg("-2").arg(&readfile[1]);
    } 

    match output.status() {
        Ok(status) if status.success() => {
            info!("SPAdes completed successfully for {:?}", binfile);
            Ok(true)
        }
        Ok(_) => {
            error!("NOTE: SPAdes failed for {:?}! This is likely due to small input bin and insufficient k-mer counts from readset."
                ,binfile.iter()
                .rev()
                .nth(2)  // nth(2) gives the third component (0-based index)
                .map(|bin| bin.to_string_lossy().to_string()).unwrap());
            Ok(false)
        }
        Err(e) => {
            error!("Error: Failed to execute SPAdes command - {}", e);
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
            Err(e) => error!("Failed to read matching file: {:?}", e),
        }
    }
    Ok(())
}

pub fn get_outputname_connected(input_file: &str) -> PathBuf {
    let path = Path::new(input_file);
    
    let output_dir = path.parent().unwrap_or_else(|| Path::new("."));
    
    let filename = path.file_stem()
        .map(|stem| stem.to_str().unwrap_or("default"))
        .unwrap_or("default");
    output_dir.join(format!("{}_connected_scaffolds", filename))
}

pub fn get_output_binname(bin_fasta: &str) -> PathBuf {
    let path = Path::new(bin_fasta);
    
    let output_dir = path.parent().unwrap_or_else(|| Path::new("."));
    let filename = path.file_stem()
        .map(|stem| stem.to_str().unwrap_or("default"))
        .unwrap_or("default");
    output_dir.join(format!("{}_enriched.fasta", filename))
}


pub fn filterscaffold(input_file: &PathBuf) -> io::Result<()> {

    let output_file = match Path::new(input_file)
        .extension() {
        Some(ext) => {
            let stem = Path::new(input_file)
                .file_stem()
                .unwrap_or_default();
            let parent = Path::new(input_file)
                .parent()
                .unwrap_or_else(||
                Path::new(""));
            parent.join(
        format!("{}_filtered.{}",
            stem.to_string_lossy(),
            ext.to_string_lossy()))
        }
        None => Path::new(input_file)
            .with_file_name(
            format!("{}_filtered",
            input_file.display())),
    };
    let input = File::open(input_file)?;
    let mut output = File::create(&output_file)?;
    let reader = BufReader::new(input);

    let mut current_header = String::new();
    let mut current_sequence = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // Process the previous sequence, if any.
            if !current_header.is_empty() && current_sequence.len() >= 500 {
                writeln!(output, "{}", current_header)?;
                writeln!(output, "{}", current_sequence)?;
            }
            current_header = line;
            current_sequence.clear();
        } else {
            current_sequence.push_str(&line.trim());
        }
    }

    if !current_header.is_empty() && current_sequence.len() >= 500 {
        writeln!(output, "{}", current_header)?;
        writeln!(output, "{}", current_sequence)?;
    }

    println!("Filtered sequences saved to: {}", output_file.to_string_lossy());

    Ok(())
}