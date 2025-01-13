use clap::{arg, Subcommand, Parser};
use std::path::PathBuf;

mod index;
mod merge;
mod gfaparser;
mod readfetch;
mod utility;
use index::indexfastqreads;
use merge::merge;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Index FASTQ files
    Index {
        /// Input directory containing FASTQ files
        #[arg(short = 'r', long = "readdir", help = "Directory containing read files")]
        readdir: PathBuf,
        /// Output directory for the database
        #[arg(short = 'o', long = "outdir", help = "Directory to output index file reads.db, Preferred to be the bin directory")]
        outdir: PathBuf,
    },
    /// Merge bins across samples
    Merge {
        /// Directory containing fasta files of bins
        #[arg(short = 'b', long = "bindir", help = "Directory containing fasta files of bins")]
        bindir: PathBuf,

        /// Directory containing gfa files for metagenomic samples (in gfa1.2 format)
        #[arg(short = 'g', long = "gfadir", help = "Directory containing gfa files")]
        gfadir: PathBuf,

        /// Directory containing sample-wise assembly contigs file in fasta format
        #[arg(short = 'a', long = "assemblydir", help = "Directory containing assembly contigs")]
        assemblydir: PathBuf,

        /// Directory containing mapids files derived from alignment sam/bam files
        #[arg(short = 'm', long = "mapdir", help = "Directory containing mapids files")]
        mapdir: PathBuf,

        /// Directory containing read files
        #[arg(short = 'r', long = "readdir", help = "Directory containing read files")]
        readdir: PathBuf,

        /// Bin file extension
        #[arg(short = 'f', long = "format", default_value = "fasta", help = "Bin file extension")]
        format: String,

        /// Number of threads to use
        #[arg(short = 't', long = "threads", default_value_t = 8, help = "Number of threads to use")]
        threads: usize,

        /// Minimum overlap length
        #[arg(short = 'l', long = "min_overlaplen", default_value_t = 1000, help = "Minimum overlap length")]
        min_overlaplen: usize,
    },
}

fn main() {
    env_logger::init();
    let cli = Cli::parse(); // Parses the command-line arguments into the `Cli` struct.
    match cli.command { // Matches the chosen subcommand.
        Commands::Index { readdir, outdir } => {
            println!("Running index command...");
            if let Err(e) = indexfastqreads(&readdir, &outdir) {
                eprintln!("Error in indexing: {}", e);
            }
        }
        Commands::Merge {
            bindir,
            gfadir,
            assemblydir,
            mapdir,
            readdir,
            format,
            threads,
            min_overlaplen,
        } => {
            println!("Running merge command...");
            if let Err(e) = merge(
                bindir,
                gfadir,
                assemblydir,
                mapdir,
                readdir,
                format,
                threads,
                min_overlaplen,
            ) {
                eprintln!("Error in merging: {}", e);
            }
        }
    }
}