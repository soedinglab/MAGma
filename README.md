# mergebins

Merge bins across samples

    Usage: mergebins [OPTIONS] --bindir <BINDIR> --gfadir <GFADIR> --assemblydir <ASSEMBLYDIR> --mapdir <MAPDIR> --readdir <READDIR>

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

