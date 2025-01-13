# mergebins

Usage:

    Index FASTQ files

    Usage: mergebins index --readdir <READDIR> --outdir <OUTDIR>

    Options:
    -r, --readdir <READDIR>  Directory containing read files
    -o, --outdir <OUTDIR>    Directory to output index file reads.db, Preferred to be the bin directory
    -h, --help               Print help

Run merge bins
    Merge bins across samples

    Usage: mergebins merge [OPTIONS] --bindir <BINDIR> --gfadir <GFADIR> --assemblydir <ASSEMBLYDIR> --mapdir <MAPDIR> --readdir <READDIR>

    Options:
    -b, --bindir <BINDIR>                  Directory containing fasta files of bins
    -g, --gfadir <GFADIR>                  Directory containing gfa files
    -a, --assemblydir <ASSEMBLYDIR>        Directory containing assembly contigs
    -m, --mapdir <MAPDIR>                  Directory containing mapids files
    -r, --readdir <READDIR>                Directory containing read files
    -f, --format <FORMAT>                  Bin file extension [default: fasta]
    -t, --threads <THREADS>                Number of threads to use [default: 8]
    -l, --min_overlaplen <MIN_OVERLAPLEN>  Minimum overlap length [default: 1000]
    -h, --help                             Print help

Note: sqlite3 should be installed