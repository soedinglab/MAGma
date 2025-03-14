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
            self.segment_to_scaffold.insert(segment.to_string(), scaffold.clone());
        }
        self.scaffold_to_segments.insert(scaffold, segments);
    }

    // Adds an overlap between two segments, linking their scaffolds
    fn add_overlap(&mut self, seg1: String, seg2: String) {
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
                    // Ignore long path for short contig, due to repeat-region or uncertain assembly structure
                    let fields: Vec<&str> = record.split('\t').collect();
                    let scaffold = fields.get(1).unwrap_or(&"");
                    let scaffold_length = scaffold.split("_length_")
                        .nth(1)
                        .and_then(|s| s.split('_').next())
                        .and_then(|s| s.parse::<usize>().ok())
                        .unwrap_or(0);

                    // Filter condition: path length should not be significantly larger than scaffold length
                    if scaffold_length >= 500 {
                        path_lines.push(record);
                    }
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
            // let segments: Vec<String> = segments
            //     .split([';', ','])
            //     .map(|s| s.to_string())
            //     .collect();
            let segments: Vec<String> = segments.split([';', ',']).map(String::from).collect();
            let first_last = vec![
                segments.first().cloned().unwrap_or_default(),
                segments.last().cloned().unwrap_or_default(),
            ];
            graph.add_path(scaffold.to_string(), first_last);
        }
    }

    // Process link lines
    for line in link_lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if let [_, seg1, dir1, seg2, dir2, _] = &fields[..] {
            graph.add_overlap(format!("{}{}", seg1, dir1), format!("{}{}", seg2, dir2));
        }
    }
    Ok(graph)
}


fn write_combined_fasta(
    sample: &String,
    enriched_scaffolds: &HashSet<String>,
    assembly_fasta: &str,
    output_fasta: PathBuf,
    create_new: bool,
) -> Result<(), Box<dyn std::error::Error>> {

    debug!("Entering write combined_fasta");
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

    debug!("Enriched scaffolds len:{}", enriched_scaffolds.len());

    let mut current_scaffold = String::new();
    let mut current_sequence = String::new();
    let mut is_first_scaffold = true;

    for line in assembly_content.lines() {
        if line.starts_with(">") {
            if !is_first_scaffold {
                if enriched_scaffolds.contains(&current_scaffold)
                    && current_sequence.len() >= 1000 {
                    writeln!(output_file, ">{}C{}", sample, current_scaffold)?;
                    writeln!(output_file, "{}", current_sequence)?;
                }
            } else {
                is_first_scaffold = false;
            }
    
            if let Some(pos) = line.find('C') {
                current_scaffold = line[pos + 1..].to_string();
            } else {
                current_scaffold = line.trim_start_matches(">").to_string();
            }
    
            current_sequence.clear();
        } else {
            current_sequence.push_str(line.trim());
        }
    }
    
    if enriched_scaffolds.contains(&current_scaffold)
        && current_sequence.len() >= 1000 {
        writeln!(output_file, ">{}C{}", sample, current_scaffold)?;
        writeln!(output_file, "{}", current_sequence)?;
    }
    Ok(())
}

pub fn read_fasta_wosid(fasta_file: &str) -> io::Result<HashSet<String>> {
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
    sample: &String,
    gfa_file: &str,
    bin_fasta: &str,
    assembly_fasta: &str,
    outputbin: PathBuf,
    create_new: bool,
) -> io::Result<HashSet<String>> {

    let graph = parse_gfa(gfa_file)?;
    debug!("obtained graph for {:?}", gfa_file);

    // eg: bin_scaffolds = {k141_1, k141_2, ..., k141_N}
    let bin_scaffolds = read_fasta_wosid(bin_fasta)?;
    debug!("bin scaffold {:?}", bin_scaffolds);

    // Get only directly linked neighbors to bin contigs
    let mut connected_scaffolds = bin_scaffolds.clone();

    for node in &bin_scaffolds {
        if let Some(neighbors) = graph.scaffold_graph.get(node) {
            connected_scaffolds.extend(neighbors.iter().cloned()); // Add neighbors to the result
        }
    }
    debug!("obtained connected components and assembly fasta is {:?}", assembly_fasta);
    debug!("output fasta {:?}", utility::get_output_binname(outputbin.to_str().expect("")));
    // eg: output_fasta = <bindir>/0_combined/S1_enriched.fasta
    let _ = write_combined_fasta(
            sample,
            &connected_scaffolds,
            assembly_fasta,
            utility::get_output_binname(outputbin.to_str().expect("")),
            create_new
        );
    // eg: output_fasta = <bindir>/bin_1_connected_components
    // let binding = utility::get_outputname_connected(outputbin.to_str().expect(""));
    // let output_file = binding
    //     .to_str()
    //     .expect("Failed to convert PathBuf to &str");

    // let mut output_cc = if create_new {
    //     File::create(&utility::get_outputname_connected(output_file))? // Create a new file if flag is true
    // } else {
    //     OpenOptions::new()
    //         .create(true)
    //         .append(true)  // Open in append mode if flag is false
    //         .open(&utility::get_outputname_connected(output_file))?
    // };

    // writeln!(output_cc, "Connected components:")?;
    // for (i, component) in components.iter().enumerate() {
    //     writeln!(output_cc, "Component {}: {:?}", i + 1, component)?;
    // }
    // debug!("connected components are written to {:?}", output_cc);
    Ok(connected_scaffolds)
}

// pub fn get_outputname_connected(input_file: &str) -> PathBuf {
//     let path = Path::new(input_file);
    
//     let output_dir = path.parent().unwrap_or_else(|| Path::new("."));
    
//     let filename = path.file_stem()
//         .map(|stem| stem.to_str().unwrap_or("default"))
//         .unwrap_or("default");
//     output_dir.join(format!("{}_connected_scaffolds", filename))
// }