# MAss

Merge and Assemble bins across samples

    Usage: mass [OPTIONS] --bindir <BINDIR> --mapdir <MAPDIR> --readdir <READDIR>

    Options:
    -b, --bindir <BINDIR>            Directory containing fasta files of bins
    -i, --ani <ANI>                  ANI for clustering bins (%) [default: 99]
    -m, --mapdir <MAPDIR>            Directory containing mapids files
    -r, --readdir <READDIR>          Directory containing read files
    -f, --format <FORMAT>            Bin file extension [default: fasta]
    -t, --threads <THREADS>          Number of threads to use [default: 8]
        --split                      Split clusters into sample-wise bins before processing
    -g, --gfadir <GFADIR>            Directory containing gfa files
    -a, --assemblydir <ASSEMBLYDIR>  Directory containing assembly contigs
        --assembler <ASSEMBLER>      assembler choice for reassembly step (spades|megahit) [default: spades]
    -h, --help                       Print help
    -V, --version                    Print version

Mapid files can be generated using aligner2counts (https://github.com/soedinglab/binning_benchmarking/tree/main/util#aligner2counts) with `only-mapids` option. 