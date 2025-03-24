# MAGma
MAGma is a tool to maximize the yield of <u>M</u>etagenome-<u>A</u>ssembled <u>G</u>enomes (MAGs) through <u>M</u>erging and re<u>A</u>ssembly.

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

Option 1: Use conda package

    conda install -c bioconda magma
    or
    mamba install -c bioconda magma // faster installation

Option 2: Use the pre-built x86-64 Linux statically compiled executable. 

    wget https://github.com/soedinglab/MAGma/releases/download/latest/magma
    chmod +x magma
    ./magma -h

To use this, [CheckM2](https://github.com/chklovski/CheckM2), [skani](https://github.com/bluenote-1577/skani), [SPAdes](https://github.com/ablab/spades) and [MEGAHIT](https://github.com/voutcn/megahit) must already be installed and available in your PATH.

Option 3: Build from source

    git clone https://github.com/soedinglab/MAGma.git
    cd MAGma
    conda env create -f environment.yml
    conda activate magma_env
    cargo install --path . 
    magma -h



### Notes
Input contigs should be header prefixed with the sample ID, separated by 'C'. Perform mapping and binning on contig files with these updated contig ids.

Mapid files can be generated using aligner2counts (https://github.com/soedinglab/binning_benchmarking/tree/main/util#aligner2counts) with `only-mapids` option.

File name: `<sampleid>_mapids`
```
read1_id    sampleidCcontig1_id
read2_id    sampleidCcontig2_id
read2_id    sampleidCcontig4_id
read3_id    sampleidCcontig2_id
read4_id    sampleidCcontig3_id
read4_id    sampleidCcontig4_id
```

If input bins are not separated by sample IDs, such as when using MetaBAT2 on a concatenated set of contigs, use the `--split` option to automatically separate clusters by sample IDs.

Make sure that headers in the read fastq files have read_id separated by space/tab (not by `.`) from other sequencer details. This is important for `seqtk` to fetch reads correctly.

`Correct header: @SRR25448374.1 A00214R:157:HLMVMDSXY:1:1101:19868:1016:N:0:CAAGTTATTG+NCGCAGAGTA.length=151#0/1`

`Doesn't work: @SRR25448374.1.A00214R:157:HLMVMDSXY:1:1101:19868:1016:N:0:CAAGTTATTG+NCGCAGAGTA.length=151#0/1`

MAGma accepts both paired-end (in separate files like SRR25448374_1.fastq and SRR25448374_2.fastq) and single-end read files.

Sample IDs must be in the file name of fastq and mapid files. (E.g., SRR25448374_1.fastq & SRR25448374_2.fastq or SRR25448374.fastq and SRR25448374_mapids)
