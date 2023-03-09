# Single cell nanoCut&Tag protocol
### Author1, Author2, Author3...
<hr>
Single cell nanoCut&Tag

Author1, Author2, Author3...

Nature Protocols 2023 (link)

# Overview
This documentation will cover the Single cell nanoCut&Tag, Nature Protocols 2023, in deepth to successfully reproduce step by step bioinformatic workflow described in the paper.
The pipeline is composed three majors axes, [Set Up](#set-up), [Preprocessing](#preprocessing) and [Downstream Analysis](#downstream-analysis) which will be broken down into sub-units to facilitate the go-through.

# Data availability
All raw and processed files can be found as supplementary files in the [GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198467).
`seurat.rds` object can also be used to start the analysis from [Downstream Analysis](#downstream-analysis).

> In this tutorail, only `GSM5949206` and `GSM5949208` will be processed

# Set up
The whole project has been run on a High Performance Computing (HPC) linux cluster under CentOS (release:7.9.2009) with htcondor workflow management system.
If you fancy using MacOS or Windows, please design your set up accordingly.

## Prepare environment
A conda environment will be used to set up an isolated architecture reducing troubleshooting.
Conda can be installed via miniconda following [miniconda guidelines](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

### Create conda environment
```
conda create --name NatProt -y -c conda-forge -c bioconda python==3.10.8 htcondor==10.2.1 snakemake==7.24.0 cookiecutter==2.1.1 git==2.39.1 sra-tools==3.0.3 r-ggplot2==3.4.1 r-argparse==2.1.5 r-funr==0.3.2 r-patchwork==1.1.2
```
Sit in your new environment till the end of the procedure.
```
conda activate NatProt
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
fasterq-dump -f -e 1 --split-files --include-technical -o SRA.fastq SRR18305888
mv SRA_1.fastq P23209_1001_S1_L001_I1_001.fastq
mv SRA_2.fastq P23209_1001_S1_L001_R1_001.fastq
mv SRA_3.fastq P23209_1001_S1_L001_R2_001.fastq
mv SRA_4.fastq P23209_1001_S1_L001_R3_001.fastq
gzip *.fastq
mkdir -p ~/NatProt/Data/fastq/P23209/P23209_1001/02-FASTQ/211025_A00187_0608_AHL2NKDRXY/
mv *.fastq.gz ~/NatProt/Data/fastq/P23209/P23209_1001/02-FASTQ/211025_A00187_0608_AHL2NKDRXY/
```
```
fasterq-dump -f -e 1 --split-files --include-technical -o SRA.fastq SRR18305889
mv SRA_1.fastq P23209_1001_S1_L002_I1_001.fastq
mv SRA_2.fastq P23209_1001_S1_L002_R1_001.fastq
mv SRA_3.fastq P23209_1001_S1_L002_R2_001.fastq
mv SRA_4.fastq P23209_1001_S1_L002_R3_001.fastq
gzip *.fastq
mkdir -p ~/NatProt/Data/fastq/P23209/P23209_1001/02-FASTQ/211025_A00187_0608_AHL2NKDRXY/
mv *.fastq.gz ~/NatProt/Data/fastq/P23209/P23209_1001/02-FASTQ/211025_A00187_0608_AHL2NKDRXY/
```
```
fasterq-dump -f -e 1 --split-files --include-technical -o SRA.fastq SRR18305884
mv SRA_1.fastq P24004_1001_S1_L001_I1_001.fastq
mv SRA_2.fastq P24004_1001_S1_L001_R1_001.fastq
mv SRA_3.fastq P24004_1001_S1_L001_R2_001.fastq
mv SRA_4.fastq P24004_1001_S1_L001_R3_001.fastq
gzip *.fastq
mkdir -p ~/NatProt/Data/fastq/P24004/P24004_1001/02-FASTQ/211221_A00621_0569_BHTNK3DRXY/
mv *.fastq.gz ~/NatProt/Data/fastq/P24004/P24004_1001/02-FASTQ/211221_A00621_0569_BHTNK3DRXY/
```
```
fasterq-dump -f -e 1 --split-files --include-technical -o SRA.fastq SRR18305885
mv SRA_1.fastq P24004_1001_S1_L002_I1_001.fastq
mv SRA_2.fastq P24004_1001_S1_L002_R1_001.fastq
mv SRA_3.fastq P24004_1001_S1_L002_R2_001.fastq
mv SRA_4.fastq P24004_1001_S1_L002_R3_001.fastq
gzip *.fastq
mkdir -p ~/NatProt/Data/fastq/P24004/P24004_1001/02-FASTQ/211221_A00621_0569_BHTNK3DRXY/
mv *.fastq.gz ~/NatProt/Data/fastq/P24004/P24004_1001/02-FASTQ/211221_A00621_0569_BHTNK3DRXY/
```

> This step takes a while...

## Clone github repository
```
cd ~/NatProt
git clone https://github.com/bartosovic-lab/single-cell-nano-cut-tag
```

> If you encounter authentication errors, you need to create a [personal access token](https://docs.github.com/fr/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)

## Install cellranger
Cellranger is a bit more complex to fit inside a conda environment and is generally heavy to store.
On most HPCs running bioinformatic pipelines, Cellranger is already installed.
However if you wish to run the pipeline on a separated workstation, you can follow [10Xgenomics guidelines](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation).

> If cellranger is not installed on your favorite HPC, please contact your HPC support.

## HPC profiles
The pipeline is implemented in workflow management software known as snakemake.
It communicates with HPCs to run paralellized jobs to speed up the process.
Like previously mentionned at the beginning of the [Set Up](#set-up), the conda environment has been built to create a communication between snakemake and htcordor scheduler, therefore, htcondor package has been installed in the conda environment.

For htcondor workflow management, we will follow this [guidelines](https://github.com/Snakemake-Profiles/htcondor)
```
mkdir -p ~/.config/snakemake
cd ~/.config/snakemake
```
```
cookiecutter https://github.com/Snakemake-Profiles/htcondor.git
```
> At `profile_name [htcondor]` press `enter` and select a path for your log files, something like `~/condor_jobs`

If your HPC is running on a different scheduler, you can install different style of profiles, like [slurm](https://github.com/Snakemake-Profiles/slurm).

> Do not forget to change the profile in your snakemake command line : 
```
--profile slurm
```
> If you do not have access to an HPC, you can simply remove the profile options from the snakemake command line.

## Modify config.yaml

```
mkdir -p ~/tmp/NatProt
```

Modify enverything by hand
For the sake of the tutorial, only 3 modalities datasets will be used bcdCT_MB21_02 and bcdCT_MB21_04

# Preprocessing
The first steps of processing from fastq files to cell picking will be done by the workflow management system, snakemake.
It will cover the following steps :
 - Modality Demultiplexing
 - Reads alignment with Cellranger
 - Transform bam files to bigwig files
 - Call peaks with Macs2
 - Label fragments with associated barcode
 - Label overlapping peaks with associated barcode
 - Output barcode metrics
 - Redo cell selection

All outputed files will be automatically generated and will be used to run the R markdown vignette below.

Run snakemake :
```
cd ~/NatProt/single-cell-nano-cut-tag
snakemake --snakefile workflow/Snakefile --cores 16 --profile htcondor -p
```

# Downstream analysis
Description of Rmarkdown

