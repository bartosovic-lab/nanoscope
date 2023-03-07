# Single cell nanoCut&Tag protocol
### Author1, Author2, Author3...
<hr>
Single cell nanoCut&Tag

Author1, Author2, Author3...

Nature Protocols 2023 (link)

# Overview
This documentation will cover the Single cell nanoCut&Tag, Nature Protocols 2023, in deepth to successfully reproduce the step by step bioinformatic workflow described in the paper.
The pipeline is composed two majors axes, [Preprocessing](#preprocessing) and [Downstream Analysis](#downstream-analysis) which will be broken down in sub-units to facilitate the go-through.

# Data Availability
All raw and processed files can be found as supplementary files in the [GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198467).
Seurat .rds object can also be used to start the analysis at the [Downstream Analysis](#downstream-analysis).

# Set Up
The whole project has been run on a High Performance Computing (HPC) linux cluster under CentOS(check version).
If you fancy using MacOS or Windows, please design your set up accordingly.

## Prepare environment
A conda environment will be used to set up an isolated architecture reducing troubleshooting.
Conda can be installed via miniconda following [miniconda guidelines](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

### Create conda environment
```
conda create --name NatProt
```
```
conda activate NatProt
```
```
conda install -y -c conda-forge -c bioconda python==3.10.8 pysam==0.20.0 htcondor==10.2.1 python-levenshtein==0.20.9 pyyaml==6.0 deeptools==3.5.1 r-base==4.2.2 snakemake==7.24.0 samtools=1.16.1 cookiecutter==2.1.1 regex==2022.10.31 gzip==1.12 contextlib2==21.6.0 bedtools==2.30.0 macs2==2.2.7.1 git==2.39.2
```

## Download the raw data
The raw data as fastq files can be downloaded throughout the [SRA-Toolkit](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump).
### Navigate to a new directory to store the data
```
mkdir -p ~/Data/NatProt
cd ~/Data/NatProt
```
### Download SRA
```
~/miniconda3/envs/NatProt/bin/fastq-dump -F --split-files SRR18305888
~/miniconda3/envs/NatProt/bin/fastq-dump -F --split-files SRR18305889
```

## Clone github repository
