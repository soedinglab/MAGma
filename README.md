# MAGma
MAGma is a tool to maximize the yield of Metagenome-Assembled Genomes (MAGs) through Merging and reAssembly.

    Usage: magma [OPTIONS] --bindir <BINDIR> --mapdir <MAPDIR> --readdir <READDIR>

    Options:
    -b, --bindir <BINDIR>
            Directory containing fasta files of bins
    -i, --ani <ANI>
            ANI for clustering bins (%) [default: 99]
    -c, --completeness <COMPLETENESS_CUTOFF>
            Minimum completeness of bins (%) [default: 50]
    -p, --purity <PURITY_CUTOFF>
            Mininum purity of bins (%) [default: 95]
    -m, --mapdir <MAPDIR>
            Directory containing mapids files
    -r, --readdir <READDIR>
            Directory containing read files
    -f, --format <FORMAT>
            Bin file extension [default: fasta]
    -t, --threads <THREADS>
            Number of threads to use [default: 8]
        --split
            Split clusters into sample-wise bins before processing
        --assembler <ASSEMBLER>
            assembler choice for reassembly step (spades|megahit) [default: spades]
    -h, --help
            Print help
    -V, --version

Mapid files can be generated using aligner2counts (https://github.com/soedinglab/binning_benchmarking/tree/main/util#aligner2counts) with `only-mapids` option.

File name: `<sampleid>_mapids`

read1_id    contig1_id
read2_id    contig2_id
read2_id    contig4_id
read3_id    contig2_id
read4_id    contig3_id
read4_id    contig4_id

If input bins are not separated by sample IDs, such as when using MetaBAT2 on a concatenated set of contigs, use the `--split` option to automatically separate clusters by sample IDs.