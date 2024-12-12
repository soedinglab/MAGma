use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::fs::{File, read_to_string, OpenOptions};
use std::io::{self, BufRead, Write};

#[derive(Debug)]
struct GfaGraph {
    // Mapping from scaffold ID to its constituent segments
    scaffold_to_segments: HashMap<String, Vec<String>>,
    // Mapping from segment ID to its scaffold
    segment_to_scaffold: HashMap<String, String>,
    // Graph of scaffolds connected by overlaps
    scaffold_graph: HashMap<String, HashSet<String>>,
}

impl GfaGraph {
    /// Creates a new empty GfaGraph.
    fn new() -> Self {
        GfaGraph {
            scaffold_to_segments: HashMap::new(),
            segment_to_scaffold: HashMap::new(),
            scaffold_graph: HashMap::new(),
        }
    }

    /// Parses a P line to map scaffolds to segments.
    fn add_path(&mut self, scaffold: String, segments: Vec<String>) {
        for segment in &segments {
            // Remove the orientation (+/-) for the segment mapping
            let clean_segment = segment.trim_end_matches(&['+', '-']);
            self.segment_to_scaffold.insert(clean_segment.to_string(), scaffold.clone());
        }
        self.scaffold_to_segments.insert(scaffold, segments);
    }

    /// Adds an overlap between two segments, linking their scaffolds.
    fn add_overlap(&mut self, seg1: String, seg2: String) {
        // Remove orientations (+/-) for clean segment names
        let seg1 = seg1.trim_end_matches(&['+', '-']).to_string();
        let seg2 = seg2.trim_end_matches(&['+', '-']).to_string();

        if let (Some(scaf1), Some(scaf2)) = (
            self.segment_to_scaffold.get(&seg1),
            self.segment_to_scaffold.get(&seg2),
        ) {
            if scaf1 != scaf2 {
                self.scaffold_graph
                    .entry(scaf1.clone())
                    .or_insert_with(HashSet::new)
                    .insert(scaf2.clone());
                self.scaffold_graph
                    .entry(scaf2.clone())
                    .or_insert_with(HashSet::new)
                    .insert(scaf1.clone());
            }
        }
    }

    /// Finds connected components of scaffolds.
    fn connected_components(&self) -> Vec<HashSet<String>> {
        let mut visited = HashSet::new();
        let mut components = Vec::new();
    
        for scaffold in self.scaffold_graph.keys() {
            if visited.contains(scaffold) {
                continue; // Skip already visited scaffolds
            }
            let mut stack = vec![scaffold.clone()];
            let mut component = HashSet::new();

            while let Some(current) = stack.pop() {
                if !visited.insert(current.clone()) {
                    continue; // Skip if already visited
                }
                component.insert(current.clone());
                if let Some(neighbors) = self.scaffold_graph.get(&current) {
                    stack.extend(neighbors.iter().filter(|n| !visited.contains(*n)).cloned());
                }
            }
            components.push(component);
        }
        components
    }
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let open = File::open(filename);
    let file = open?;
    Ok(io::BufReader::new(file).lines())
}

/// Parses a GFA file and constructs a GfaGraph.
fn parse_gfa(file_path: &str) -> io::Result<GfaGraph> {

    let mut graph = GfaGraph::new();
    let mut path_lines = Vec::new();
    let mut link_lines = Vec::new();

    if let Ok(lines) = read_lines(file_path) {
        for line in lines {
            if let Ok(record) = line {
                if record.starts_with('P') {
                    path_lines.push(record);
                } else if record.starts_with('L') {
                    link_lines.push(record);
                }
            }
        }
    }

    // Step 2: Process Path lines first
    for line in path_lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if let [_, scaffold, segments, _] = &fields[..] {
            let segments: Vec<String> = segments
                .split([';', ','])
                .map(|s| s.to_string())
                .collect();
            graph.add_path(scaffold.to_string(), segments);
        }
    }

    for line in link_lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if let [_, seg1, _, seg2, _, _] = &fields[..] {
            graph.add_overlap(seg1.to_string(), seg2.to_string());
        }
    }
    Ok(graph)
}


fn get_output_filename(input_file: &str) -> PathBuf {
    let path = Path::new(input_file);
    
    // Get the directory of the GFA file (unwrap safely with a default)
    let output_dir = path.parent().unwrap_or_else(|| Path::new("."));
    
    let filename = path.file_stem()
        .map(|stem| stem.to_str().unwrap_or("default"))
        .unwrap_or("default");
    output_dir.join(format!("{}_connected_scaffolds", filename))
}

fn get_output_scaffoldname(bin_fasta: &str) -> PathBuf {
    let path = Path::new(bin_fasta);
    
    // Get the directory of the GFA file (unwrap safely with a default)
    let output_dir = path.parent().unwrap_or_else(|| Path::new("."));
    
    let filename = path.file_stem()
        .map(|stem| stem.to_str().unwrap_or("default"))
        .unwrap_or("default");
    output_dir.join(format!("{}_enriched.fasta", filename))
}

fn write_combined_fasta(
    bin_scaffolds: &HashSet<String>,
    connected_scaffolds: &HashSet<String>,
    assembly_fasta: &str,
    output_fasta: PathBuf,
    create_new: bool,
) -> io::Result<HashSet<String>> {
    let assembly_content = read_to_string(assembly_fasta)?;
    let output_file = if create_new {
        File::create(&output_fasta)?
    } else {
        OpenOptions::new()
            .create(true)
            .append(true)
            .open(&output_fasta)?
    };

    let mut output_file = io::BufWriter::new(output_file);

    let mut enriched_scaffolds = bin_scaffolds.clone();
    // a complete set of scaffolds
    enriched_scaffolds.extend(connected_scaffolds.iter().cloned());

    let mut current_scaffold = String::new();
    let mut current_sequence = String::new();

    for line in assembly_content.lines() {
        if line.starts_with(">") {
            if !current_scaffold.is_empty() 
            && enriched_scaffolds.contains(&current_scaffold)
            && current_sequence.len() >=300 {
                writeln!(output_file, ">{}", current_scaffold)?;
                writeln!(output_file, "{}", current_sequence)?;
            }
            current_scaffold = line.trim_start_matches(">").to_string();
            current_sequence.clear();
        } else {
            current_sequence.push_str(line);
        }
    }

    // Last line scaffold
    if !current_scaffold.is_empty() 
        && enriched_scaffolds.contains(&current_scaffold)
        && current_sequence.len() >=300 {
        writeln!(output_file, ">{}", current_scaffold)?;
        writeln!(output_file, "{}", current_sequence)?;
    }

    Ok(enriched_scaffolds)
}

fn read_fasta(fasta_file: &str) -> io::Result<HashSet<String>> {
    let content = read_to_string(fasta_file)?;
    let mut scaffolds = HashSet::new();
    for line in content.lines() {
        if line.starts_with(">") {
            let scaffold_name = line.trim_start_matches(">")
            .split_whitespace()
            .next()
            .unwrap_or(line.trim_start_matches(">"))
            .split_once('C') // Split into two parts at the first 'C'
            .map(|(_, after_c)| after_c) // Take the part after 'C'
            .unwrap_or("");
            scaffolds.insert(scaffold_name.to_string());
        }
    }
    Ok(scaffolds)
}

fn read_mapid_file(mapid_file: &str) -> io::Result<HashMap<String, String>> {
    let mut read_to_scaffold = HashMap::new();
    let file = File::open(mapid_file)?;
    let reader = io::BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 2 {
            let read_id = parts[0].to_string();
            let scaffold_id = parts[1].to_string();
            read_to_scaffold.insert(read_id, scaffold_id);
        }
    }
    Ok(read_to_scaffold)
}

/// Writes selected reads to a FASTQ file.
fn write_selected_reads(
    fastq_file: &str,
    enriched_scaffolds: &HashSet<String>,
    mapid_file: &str,
    output_fastq: PathBuf,
    create_new: bool,
) -> io::Result<()> {
    let read_to_scaffold = read_mapid_file(mapid_file)?;
    let selected_reads: HashSet<String> = read_to_scaffold
        .iter()
        .filter(|(_, scaffold)| enriched_scaffolds.contains(*scaffold))
        .map(|(read, _)| read.clone())
        .collect();
    let file = File::open(fastq_file)?;
    let reader = io::BufReader::new(file);
    let output_file = if create_new {
        File::create(&output_fastq)? // Create a new file if flag is true
    } else {
        OpenOptions::new()
            .create(true)
            .append(true)  // Open in append mode if flag is false
            .open(&output_fastq)?
    };

    let mut output_file = io::BufWriter::new(output_file);

    let mut lines = reader.lines();
    while let Some(header) = lines.next() {
        let sequence = lines.next().unwrap_or(Ok(String::new()))?;
        let plus_line = lines.next().unwrap_or(Ok(String::new()))?;
        let quality = lines.next().unwrap_or(Ok(String::new()))?;
        let read_id = header?.to_string();
        let read_id_trimmed = read_id
            .trim_start_matches('@')
            .trim_end_matches("/1")
            .trim_end_matches("/2")
            .to_string();
        if selected_reads.contains(&read_id_trimmed) {
            writeln!(output_file, "{}", read_id)?;
            writeln!(output_file, "{}", sequence)?;
            writeln!(output_file, "{}", plus_line)?;
            writeln!(output_file, "{}", quality)?;
        }
    }

    Ok(())
}


/// Pub function to process a GFA file and print connected components.
pub fn parse_gfa_fastq(
    gfa_file: &str,
    bin_fasta: &str,
    assembly_fasta: &str,
    mapids: &str,
    read_fastq: &str,
    outputbin: PathBuf,
    create_new: bool,
) -> io::Result<()> {

    let graph = parse_gfa(gfa_file)?;

    let bin_scaffolds = read_fasta(bin_fasta)?;
    // Find connected components of scaffolds
    let components = graph.connected_components();

    let mut connected_scaffolds = HashSet::new();
    for component in &components {
        for scaffold in component {
            if bin_scaffolds.contains(scaffold) {
                connected_scaffolds.extend(component.iter().cloned());
                break;
            }
        }
    }
    let enriched_scaffolds = write_combined_fasta(
        &bin_scaffolds, &connected_scaffolds, 
        assembly_fasta, 
        get_output_scaffoldname(outputbin.to_str().expect("")),
        create_new)?;

    let output_fastq = PathBuf::from(format!("{}",
        get_output_scaffoldname(outputbin.to_str().expect(""))
        .to_str()
        .expect("Invalid UTF-8 in file path")
        .replace(".fasta", ".fastq")));

    let _ = write_selected_reads(
        &read_fastq,
        &enriched_scaffolds,
        mapids,
        output_fastq,
        create_new);

    let binding = get_output_filename(outputbin.to_str().expect(""));
    let output_file = binding
        .to_str()
        .expect("Failed to convert PathBuf to &str");

    let mut output_cc = if create_new {
        File::create(&get_output_filename(output_file))? // Create a new file if flag is true
    } else {
        OpenOptions::new()
            .create(true)
            .append(true)  // Open in append mode if flag is false
            .open(&get_output_filename(output_file))?
    };

    writeln!(output_cc, "Connected components:")?;
    for (i, component) in components.iter().enumerate() {
        writeln!(output_cc, "Component {}: {:?}", i + 1, component)?;
    }
    println!("Connected components written to {:?}",output_file);
    Ok(())
}
