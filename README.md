# Single cell nanoCut&Tag protocol
### Author1, Author2, Author3...
<hr>
Single cell nanoCut&Tag

Author1, Author2, Author3...

Nature Protocols 2023 (link)

# Overview
This documentation will cover the Single cell nanoCut&Tag, Nature Protocols 2023, in deepth to successfully reproduce the step by step bioinformatic workflow described in the paper.
The pipeline is composed three majors axes, [Set Up](#set-up), [Preprocessing](#preprocessing) and [Downstream Analysis](#downstream-analysis) which will be broken down in sub-units to facilitate the go-through.

# Data availability
All raw and processed files can be found as supplementary files in the [GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198467).
Seurat .rds object can also be used to start the analysis at the [Downstream Analysis](#downstream-analysis).

# Set up
The whole project has been run on a High Performance Computing (HPC) linux cluster under CentOS (release:7.9.2009) with htcondor workflow management system.
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
mkdir -p ~/NatProt/Data
cd ~/NatProt/Data
```
### Download SRA
```
~/miniconda3/envs/NatProt/bin/fastq-dump -F --split-files SRR18305888
~/miniconda3/envs/NatProt/bin/fastq-dump -F --split-files SRR18305889
```

## Clone github repository
```
cd ~/NatProt
~/miniconda3/envs/NatProt/bin/git clone https://github.com/bartosovic-lab/single-cell-nano-cut-tag
```

## Install cellranger
Cellranger is a bit more complex to fit inside a conda environment and is generally heavy to store.
On most of the HPCs running bioinformatic pipelines, Cellranger is already installed.
However if you wish to run the pipeline on a separated workstation, you can follow [10Xgenomics guidelines](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation).

> If cellranger is not installed on your favorite HPC, please contact your HPC support.

## HPC profiles
The pipeline is implemented in workflow management software known as snakemake.
It communicates with HPCs to run paralellized jobs to speed up the process.
Like previously mentionned at the beginning of the [Set Up](#set-up), the conda environment has been built to create a communication between snakemake and htcordor scheduler, therefore, htcondor package has been installed in the conda environment.

If your HPC is running on a different scheduler, you can install different style of profiles, like [slurm](https://github.com/Snakemake-Profiles/slurm).

> Do not forget to change the profile in your snakemake command line : 
```
--profile slurm
```
> If you do not have access to an HPC, you can simply remove the profile options from the snakemake command line.

## Modify config.yaml

Bash script with user input or hardcoded in the config.yaml file ?

# Preprocessing
From fastq to Seurat objects for Downstream analysis
## Demultiplexing
~/miniconda3/envs/NatProt/bin/snakemake --snakefile ~/single-cell-nano-cut-tag/workflow/Snakefile_demultiplexing.smk --cores 16 --profile htcondor -p
## Cellranger
~/miniconda3/envs/NatProt/bin/snakemake --snakefile ~/single-cell-nano-cut-tag/workflow/Snakefile_cellranger.smk --cores 16 --profile htcondor -p
## Peaks calling
~/miniconda3/envs/NatProt/bin/snakemake --snakefile ~/single-cell-nano-cut-tag/workflow/Snakefile_peaks_calling.smk --cores 16 --profile htcondor -p
## Cell picking
~/miniconda3/envs/NatProt/bin/snakemake --snakefile ~/single-cell-nano-cut-tag/workflow/Snakefile_cell_picking.smk --cores 16 --profile htcondor -p
## Seurat objects
~/miniconda3/envs/NatProt/bin/snakemake --snakefile ~/single-cell-nano-cut-tag/workflow/Snakefile_seurat_objects.smk --cores 16 --profile htcondor -p
## Merge fragments
~/miniconda3/envs/NatProt/bin/snakemake --snakefile ~/single-cell-nano-cut-tag/workflow/Snakefile_merge_fragments.smk --cores 16 --profile htcondor -p
## Peak calling merge
~/miniconda3/envs/NatProt/bin/snakemake --snakefile ~/single-cell-nano-cut-tag/workflow/Snakefile_peaks_calling_merged.smk --cores 16 --profile htcondor -p
## Seurat object merged
~/miniconda3/envs/NatProt/bin/snakemake --snakefile ~/single-cell-nano-cut-tag/workflow/Snakefile_seurat_objects_merged.smk --cores 16 --profile htcondor -p

# Downstream analysis
Description of Rmarkdown

