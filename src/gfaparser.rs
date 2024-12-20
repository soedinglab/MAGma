use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::fs::{File, read_to_string, OpenOptions};
use std::io::{self, BufRead, Write};
use log::debug;
use crate::utility;

#[derive(Debug)]
struct GfaGraph {
    scaffold_to_segments: HashMap<String, Vec<String>>,
    segment_to_scaffold: HashMap<String, String>,
    scaffold_graph: HashMap<String, HashSet<String>>,
}

impl GfaGraph {
    fn new() -> Self {
        GfaGraph {
            scaffold_to_segments: HashMap::new(),
            segment_to_scaffold: HashMap::new(),
            scaffold_graph: HashMap::new(),
        }
    }

    // Parsing Path lines to map scaffolds to segments
    fn add_path(&mut self, scaffold: String, segments: Vec<String>) {
        for segment in &segments {
            let clean_segment = segment.trim_end_matches(&['+', '-']);
            self.segment_to_scaffold.insert(clean_segment.to_string(), scaffold.clone());
        }
        self.scaffold_to_segments.insert(scaffold, segments);
    }

    // Adds an overlap between two segments, linking their scaffolds
    fn add_overlap(&mut self, seg1: String, seg2: String) {
        let seg1 = seg1.trim_end_matches(&['+', '-']).to_string();
        let seg2 = seg2.trim_end_matches(&['+', '-']).to_string();
        debug!("seg1 {} and seg2 {}", seg1, seg2);
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

    // Finds connected components of scaffolds
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
    let open = File::open(&filename).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to open file {:?}: {}", filename.as_ref(), e),
        )
    })?;
    
    Ok(io::BufReader::new(open).lines())
}

fn parse_gfa(file_path: &str) -> io::Result<GfaGraph> {

    let mut graph = GfaGraph::new();
    let mut path_lines = Vec::new();
    let mut link_lines = Vec::new();

    if let Err(e) = read_lines(file_path).and_then(|lines| {
        for line in lines {
            if let Ok(record) = line {
                if record.starts_with('P') {
                    path_lines.push(record);
                } else if record.starts_with('L') {
                    link_lines.push(record);
                }
            }
        }
        Ok(())
    }) {
        eprintln!("Failed to process file {}: {}", file_path, e);
    }

    // Process Path lines first
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

    // Process link lines
    for line in link_lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if let [_, seg1, _, seg2, _, _] = &fields[..] {
            graph.add_overlap(seg1.to_string(), seg2.to_string());
            debug!("overlapping segment seg1 {} seg2 {} ", seg1.to_string(), seg2.to_string());
        }
    }
    Ok(graph)
}


fn write_combined_fasta(
    bin_scaffolds: &HashSet<String>,
    connected_scaffolds: &HashSet<String>,
    assembly_fasta: &str,
    output_fasta: PathBuf,
    create_new: bool,
) -> io::Result<HashSet<String>> {
    let assembly_content = read_to_string(assembly_fasta)?;
    debug!("read assembly content");
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

    enriched_scaffolds.extend(connected_scaffolds.iter().cloned());

    let mut current_scaffold = String::new();
    let mut current_sequence = String::new();

    for line in assembly_content.lines() {
        if line.starts_with(">") {
            debug!("assembly file header {}", line);
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

// TODO: get enriched scaffolds from the source samples but reads obtain from all samples
pub fn parse_gfa_fastq(
    gfa_file: &str,
    bin_fasta: &str,
    assembly_fasta: &str,
    outputbin: PathBuf,
    create_new: bool,
) -> io::Result<HashSet<String>> {

    let graph = parse_gfa(gfa_file)?;
    debug!("obtained graph for {:?}", gfa_file);

    // eg: bin_scaffolds = {S1Ck141_1, S2Ck141_2 ... S2Ck141_N}
    let bin_scaffolds = read_fasta(bin_fasta)?;
    debug!("bin scaffold {:?}", bin_scaffolds);

    // Find connected components of scaffolds
    // eg: components = Vec<{k141_1, k141_4}, {k141_6, k141_12, k141_564},...>
    let components = graph.connected_components();

    let mut connected_scaffolds = HashSet::new();
    for component in &components {
        for scaffold in component {
            debug!("scaffold {}", scaffold);
            if bin_scaffolds.contains(scaffold) {
                connected_scaffolds.extend(component.iter().cloned());
                break;
            }
        }
    }

    debug!("obtained connected components and assembly fasta is {:?}", assembly_fasta);
    debug!("output fasta {:?}", utility::get_output_scaffoldname(outputbin.to_str().expect("")));
    // eg: output_fasta = <bindir>/bin_1/S1_enriched.fasta
    let enriched_scaffolds = 
        write_combined_fasta(
            &bin_scaffolds,
            &connected_scaffolds, 
            assembly_fasta,
            utility::get_output_scaffoldname(outputbin.to_str().expect("")),
            create_new
        )?;

    // eg: output_fasta = <bindir>/bin_1_connected_components
    let binding = utility::get_output_filename(outputbin.to_str().expect(""));
    let output_file = binding
        .to_str()
        .expect("Failed to convert PathBuf to &str");

    let mut output_cc = if create_new {
        File::create(&utility::get_output_filename(output_file))? // Create a new file if flag is true
    } else {
        OpenOptions::new()
            .create(true)
            .append(true)  // Open in append mode if flag is false
            .open(&utility::get_output_filename(output_file))?
    };

    writeln!(output_cc, "Connected components:")?;
    for (i, component) in components.iter().enumerate() {
        writeln!(output_cc, "Component {}: {:?}", i + 1, component)?;
    }
    debug!("connected components are written to {:?}", output_cc);
    Ok(enriched_scaffolds)
}
