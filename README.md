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


# Example

    magma -b binsdir -m mapid_dir -r readdir -f fasta -t 24
    magma -b binsdir -m mapid_dir -r readdir -f fasta -t 24 --split // if input bins are not already split by sample id 



# Install
### Prerequisites

- **Rust**: Follow the instructions [here](https://www.rust-lang.org/tools/install) to install Rust.
- **Conda**: You can install Conda via [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution).


Option 2: Build from source
    git clone https://github.com/soedinglab/MAGma.git
    cd MAGma
    conda env create -f environment.yml
    conda activate magma_env
    cargo install --path . 
    magma -h


### Notes
Mapid files can be generated using aligner2counts (https://github.com/soedinglab/binning_benchmarking/tree/main/util#aligner2counts) with `only-mapids` option.

File name: `<sampleid>_mapids`
```
read1_id    contig1_id
read2_id    contig2_id
read2_id    contig4_id
read3_id    contig2_id
read4_id    contig3_id
read4_id    contig4_id
```

If input bins are not separated by sample IDs, such as when using MetaBAT2 on a concatenated set of contigs, use the `--split` option to automatically separate clusters by sample IDs.

Make sure that headers in the read fastq files have read_id separated by space/tab (not by `.`) from other sequencer details. This is important for `seqtk` to fetch reads correctly.

`Correct header: @SRR25448374.1 A00214R:157:HLMVMDSXY:1:1101:19868:1016:N:0:CAAGTTATTG+NCGCAGAGTA.length=151#0/1`

`Doesn't work: @SRR25448374.1.A00214R:157:HLMVMDSXY:1:1101:19868:1016:N:0:CAAGTTATTG+NCGCAGAGTA.length=151#0/1`

fastq and mapid files must have sampleid in the file name. (E.g., SRR25448374_1.fastq & SRR25448374_2.fastq or SRR25448374.fastq and SRR25448374_mapids)
