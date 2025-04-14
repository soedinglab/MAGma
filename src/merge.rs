use bio::io::fasta;
use std::io::{BufWriter, Write};
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::fs::{remove_file, File};
use std::io::{self, BufRead, BufReader};
use petgraph::graph::Graph;
use petgraph::graph::NodeIndex;
use petgraph::visit::Dfs;
use petgraph::{Undirected, prelude};
use std::process::Command as ProcessCommand;
use log::{debug, error, info, warn};
use glob::glob;

use crate::assess::BinQuality;

/// Compute all-vs-all ANI among bins
pub fn calc_ani(
    bins: &PathBuf,
    bin_qualities: &HashMap<String, BinQuality>,
    format: &String,
    ani_cutoff: f64,
    contamination_cutoff: f64
) -> Result<(Graph<String, (), petgraph::Undirected>, HashMap<(String, String), f64>), io::Error> {
    let ani_output: PathBuf = bins.join("ani_edges");
    let bin_files: Vec<String> = glob(&format!("{}/*.{}", bins.display(), format))
    .expect("Failed to read glob pattern")
    .filter_map(Result::ok)
    .map(|path| path.to_string_lossy().into_owned())
    .collect();

    if bin_files.is_empty() {
        error!("No fasta files found in {:?}", bins);
        return Err(io::Error::new(io::ErrorKind::NotFound, "No fasta files found"));
    }
    
    let _ = get_ani(bin_files, &ani_output);

    let mut bin_name_to_node: HashMap<String, prelude::NodeIndex> = HashMap::new();
    let file: File = File::open(ani_output.clone())?;
    let reader: BufReader<File> = io::BufReader::new(file);

    let mut graph: Graph<String, (), Undirected> = Graph::default();

    for (bin, q) in bin_qualities {
        if q.contamination < contamination_cutoff && q.completeness > 20.0 {
            bin_name_to_node
                .entry(bin.clone())
                .or_insert_with(|| graph.add_node(bin.clone()));
        }
    }

    let mut ani_details = HashMap::<(String, String), f64>::new();

    // Create a graph by add edge when ANI > ANI_threshold
    // When file is empty, no edge is formed and all nodes will be Singleton clusters.   
    for line in reader.lines().skip(1) {
        let line = line?;
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 3 {
            continue;
        }
    
        let bin1 = Path::new(columns[0])
            .file_stem()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| columns[0].to_string());
    
        let bin2 = Path::new(columns[1])
            .file_stem()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| columns[1].to_string());
    
        let ani: f64 = columns[2].parse().unwrap_or(0.0);
    
        // Skani reports pairs only if ANI is >= 80%
        ani_details.insert((bin1.clone(), bin2.clone()), ani);
        if ani < ani_cutoff {
            continue;
        }
    
        if let (Some(&node1), Some(&node2)) = (bin_name_to_node.get(&bin1), bin_name_to_node.get(&bin2)) {
            graph.add_edge(node1, node2, ());
        }
    }

    // Remove skani output file
    // remove_file(&ani_output).ok();

   Ok((graph, ani_details))
}

/// Find single-linkage connected components
pub fn get_connected_samples(
    graph: &Graph<String, (), Undirected>,
    ani_details: &HashMap<(String, String), f64>,
    ani_cutoff: f64,
) -> Vec<HashSet<String>> {
    let mut visited = HashSet::new();
    let mut connected_components = Vec::new();

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
            connected_components.push(component);
        }
    }
    let mut connected_samples: Vec<HashSet<String>> = vec![];
    for component in connected_components {
        if component.len() <=2 {
            connected_samples.push(component);
        } else {
            let mut subclusters = split_component_into_cliques(component, ani_details, ani_cutoff);
            connected_samples.append(&mut subclusters);
        }
    }
    connected_samples
}

// Get clique clusters
fn split_component_into_cliques(
    component: HashSet<String>,
    ani_details: &HashMap<(String, String), f64>,
    ani_cutoff: f64,
) -> Vec<HashSet<String>> {
    let mut subgraph = Graph::<String, (), petgraph::Undirected>::default();
    let mut node_map = HashMap::new();

    for bin in &component {
        let node_index = subgraph.add_node(bin.clone());
        node_map.insert(bin.clone(), node_index);
    }

    // Add edges if ANI is above cutoff
    for (bin1, bin2) in component.iter().flat_map(|b1| {
        component.iter().map(move |b2| (b1.clone(), b2.clone()))
    }) {
        if bin1 != bin2 {
            if let Some(&ani) = ani_details
                    .get(&(bin1.clone(), bin2.clone()))
                    .or_else(|| ani_details.get(&(bin2.clone(), bin1.clone()))) {
                if ani >= ani_cutoff {
                    if let (Some(&n1), Some(&n2)) = (node_map.get(&bin1), node_map.get(&bin2)) {
                        subgraph.add_edge(n1, n2, ());
                    }
                }
            }
        }
    }

    let mut cliques = Vec::new();
    let all_nodes: HashSet<NodeIndex> = subgraph.node_indices().collect();

    // Get maximal cliques
    bron_kerbosch(&subgraph, &mut HashSet::new(), &mut all_nodes.clone(), &mut HashSet::new(), &mut cliques);

    let mut subclusters: Vec<HashSet<String>> = cliques
        .into_iter()
        .map(|clique| clique.into_iter().collect())
        .collect();

    if subclusters.is_empty() {
        debug!("No clique was found for component {:?}. How can this happen? investigate", component);
        return component.into_iter().map(|node| vec![node].into_iter().collect()).collect();
    }

    let assigned_nodes: HashSet<String> = subclusters.iter().flatten().cloned().collect();
    let unassigned_nodes: Vec<String> = component.difference(&assigned_nodes).cloned().collect();
    let mut singletons = Vec::new();
    for node in unassigned_nodes {
        let mut potential_clique: Option<usize> = None;
        let mut valid = true;

        for (i, clique) in subclusters.iter().enumerate() {
            let mut all_ani_above_99 = true;
            let mut has_any_link = false;

            for clique_node in clique {
                if let Some(&ani) = ani_details
                        .get(&(node.clone(), clique_node.clone()))
                        .or_else(|| ani_details.get(&(clique_node.clone(), node.clone()))) {
                    if ani < 99.0 {
                        all_ani_above_99 = false;
                    } else {
                        has_any_link = true;
                    }
                }
            }

            // The node should only have links above 99% to this clique, and not to multiple cliques
            if has_any_link && all_ani_above_99 {
                if potential_clique.is_some() {
                    valid = false;
                    debug!("Node {} connects to multiple cliques", node);
                    break;
                }
                potential_clique = Some(i);
            }
        }

        if valid {
            if let Some(i) = potential_clique {
                subclusters[i].insert(node);
            } else {
                singletons.push(vec![node].into_iter().collect());
            }
        } else {
            singletons.push(vec![node].into_iter().collect());
        }
    }
    subclusters.extend(singletons);

    subclusters
}

// Detect maximal cliques
fn bron_kerbosch(
    graph: &Graph<String, (), petgraph::Undirected>,
    r: &mut HashSet<NodeIndex>,
    p: &mut HashSet<NodeIndex>,
    x: &mut HashSet<NodeIndex>,
    cliques: &mut Vec<HashSet<String>>,
) {
    if p.is_empty() && x.is_empty() {
        cliques.push(r.iter().map(|&n| graph[n].clone()).collect());
        return;
    }

    let pivot = *p.iter().chain(x.iter()).next().unwrap(); // Choose a pivot
    let neighbors: HashSet<NodeIndex> = graph
        .neighbors(pivot)
        .collect();

    let candidates: Vec<NodeIndex> = p.difference(&neighbors).copied().collect();
    
    for v in candidates {
        let mut new_r = r.clone();
        let mut new_p = p.clone();
        let mut new_x = x.clone();

        new_r.insert(v);
        new_p.retain(|&n| graph.contains_edge(v, n));
        new_x.retain(|&n| graph.contains_edge(v, n));

        bron_kerbosch(graph, &mut new_r, &mut new_p, &mut new_x, cliques);

        p.remove(&v);
        x.insert(v);
    }
}

/// Merged bin files within the cluster
pub fn combine_fastabins(
    inputdir: &Path,
    bin_samplenames: &HashSet<String>,
    combined_bins: &Path,
    format: &String,
) -> io::Result<()> {
    // Combine bins fasta into a single file
    let mut output_writer = File::create(
        combined_bins
        .join("combined.fasta"))?;
    for bin_samplename in bin_samplenames {
        let bin_file_path = inputdir.join(format!("{}.{}", bin_samplename, format));
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

/// Dereplicate final bins to remove any redundant bins
pub fn drep_finalbins(
    result_dir: &PathBuf,
    bin_qualities: &HashMap<String, BinQuality>,
    ani_cutoff: f64
) -> Result<(), Box<dyn std::error::Error>> {
    let ani_output: PathBuf = result_dir.join("ani_edges");
    let finalbin_files: Vec<PathBuf> = glob(&format!("{}/*.fasta", result_dir.display()))
        .expect("Failed to read glob pattern")
        .filter_map(Result::ok)
        .collect();

    let _ = get_ani(
        finalbin_files.iter().map(|p| p.to_string_lossy().into_owned()).collect(), 
        &ani_output
    );

    debug!("Bin qualities length after reassembly: {}", bin_qualities.len());

    let file: File = File::open(ani_output.clone())?;
    let reader: BufReader<File> = io::BufReader::new(file);

    let mut bins_to_remove: HashSet<String> = HashSet::new();
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
    
        let ani: f64 = columns[2].parse().expect("Failed to parse ANI value as float from column 3");

        if ani >= ani_cutoff {
            if let (Some(q1), Some(q2)) = (bin_qualities.get(&bin1), bin_qualities.get(&bin2)) {
                let score1 = q1.completeness - (5.0 * q1.contamination);
                let score2 = q2.completeness - (5.0 * q2.contamination);
    
                let worse_bin = if score1 > score2 {
                    &bin2
                } else if score1 < score2 {
                    &bin1
                } else if q1.contamination < q2.contamination {
                    &bin2
                } else {
                    &bin1
                };
                bins_to_remove.insert(worse_bin.clone());
            }
        }
    }

    debug!("Length of list with bins to remove: {}", bins_to_remove.len());
    let bin_names: HashSet<String> = finalbin_files
        .iter()
        .filter_map(|file| file.file_stem()?.to_str().map(|s| s.to_string()))
        .collect();

    let filtered_bin_names: HashSet<String> = bin_names
        .difference(&bins_to_remove)
        .cloned()
        .collect();

    // Remove redundant bins
    for bin in &bins_to_remove {
        let bin_file_path = result_dir.join(format!("{}.fasta", bin));
    
        if bin_file_path.exists() {
            remove_file(&bin_file_path).ok();
        }
    }
    
    if !cfg!(debug_assertions) {
        if let Err(e) = remove_file(&ani_output) {
            warn!("Failed to delete folder {:?}: {}", ani_output, e);
        }
    }
    
    // Write quality measures of bins
    let output_file_path = result_dir.join("bins_checkm2_qualities.tsv");
    let output_file = File::create(&output_file_path)?;
    let mut writer = BufWriter::new(output_file);

    writeln!(writer, "#Bin\tCompleteness\tContamination")?;
    let mut buffer = String::with_capacity(1024 * 1024);
    for (bin, quality) in bin_qualities.iter() {
        if filtered_bin_names.contains(bin) {
            buffer.push_str(&format!("{}\t{}\t{}\n", bin, quality.completeness, quality.contamination));
        }
    }
    writer.write_all(buffer.as_bytes())?;
    info!("Quality values of bins are written to {:?}", output_file_path);
    Ok(())
}


// Run skani
fn get_ani (
    inputbins:Vec<String>,
    ani_output: &PathBuf
) -> Result<(), io::Error> {
    let mut command = ProcessCommand::new("skani");
    command.arg("triangle");
    command.args(&inputbins);
    command.arg("-E");
    
    let output = command.output()?;
    
    if !output.status.success() {
        return Err(io::Error::new(io::ErrorKind::Other, "skani triangle failed"));
    }
    std::fs::write(ani_output, output.stdout)?;
    Ok(())
}
