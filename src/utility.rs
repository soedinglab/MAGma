use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::fs::{self, read_to_string, File};
use std::io::{self, BufRead, BufReader, Write};
use std::process::exit;
use log::{debug, error};
use glob::glob;

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

pub fn get_sample_names(bindir: &Path, extension: &str) -> io::Result<HashMap<String, String>> {
    let mut bin_sample_map = HashMap::new();
    for entry in fs::read_dir(bindir)? {
        let path = entry?.path();
        if !path.is_file() || path.extension().and_then(|e| e.to_str()) != Some(extension) {
            continue;
        }
        // Open the file and read the first line
        let file_name = path.file_stem()
            .and_then(|name| name.to_str())  // Convert OsStr to &str
            .unwrap_or("unknown")            // Fallback if conversion fails
            .to_string();
        let file = fs::File::open(&path)?;
        let mut reader = io::BufReader::new(file);
        let mut first_line = String::new();

        if reader.read_line(&mut first_line)? == 0 {
            continue;
        }

        let sample_id = first_line.trim()
            .trim_start_matches('>')
            .split('C')
            .next()
            .map(|s| s.to_string())
            .unwrap_or_default();

        bin_sample_map.insert(file_name, sample_id);
    }
    Ok(bin_sample_map)
}


pub fn splitbysampleid(
    bin: &PathBuf,
    bin_name: &str,
    binspecificdir: &PathBuf,
    format: &String,
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
            ensure_writer(&current_sample_id, bin_name, binspecificdir, format, &mut writers)?;
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
    format: &String,
    writers: &mut HashMap<String, File>,
) -> io::Result<()> {
    if !writers.contains_key(sample_id) {
        let output_filename = 
            binspecificdir
            .join(
            format!("{}_{}.{}"
            ,bin_name,sample_id, format));
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

pub fn read_fasta(fasta_file: &str) -> io::Result<HashSet<String>> {
    let content = read_to_string(fasta_file)?;
    let mut scaffolds = HashSet::new();
    for line in content.lines() {
        if line.starts_with(">") {
            let scaffold_name = line.trim_start_matches(">")
                .split_whitespace()
                .next()
                .unwrap_or(line.trim_start_matches(">"));
            scaffolds.insert(scaffold_name.to_string());
        }
    }
    Ok(scaffolds)
}

pub fn get_output_binname(bin_fasta: &str) -> PathBuf {
    let path = Path::new(bin_fasta);
    
    let output_dir = path.parent().unwrap_or_else(|| Path::new("."));
    let filename = path.file_stem()
        .map(|stem| stem.to_str().unwrap_or("default"))
        .unwrap_or("default");
    output_dir.join(format!("{}.fasta", filename))
}

// pub fn remove_files_matching_pattern(
//     directory: &Path,
//     pattern: &str
// ) -> io::Result<()> {
//     let search_pattern = 
//         directory
//         .join(pattern)
//         .to_string_lossy()
//         .to_string();
//     for entry in glob(&search_pattern)
//         .map_err(|e
//         | io::Error::new(
//         io::ErrorKind::Other,
//         e.to_string()))? {
//         match entry {
//             Ok(path) => {
//                 let _ = fs::remove_file(&path).is_err();
//             }
//             Err(e) => error!("Failed to read matching file: {:?}", e),
//         }
//     }
//     Ok(())
// }
