# Nanoscope 
# Analysis pipeline for single-cell nano-CUT&Tag data analysis
### Federico Ansaloni, Bastien Herve, Marek Bartosovic
<hr>

# Overview
This documentation will cover the Single cell nanoCut&Tag, Nature Protocols 2023, in deepth to successfully reproduce step by step bioinformatic workflow described in the paper.
The pipeline is composed of three majors axes, [Set Up](#set-up), [Preprocessing](#preprocessing) and [Downstream Analysis](#downstream-analysis) which will be broken down into sub-units to facilitate the go-through.

# Data availability
All raw and processed files can be found as supplementary files in the [GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198467).
`seurat.rds` object can also be used to start the analysis from [Downstream Analysis](#downstream-analysis).

> In this tutorial, only `GSM5949206` and `GSM5949208` will be processed

# Set up
The whole project has been run on a High Performance Computing (HPC) linux cluster under CentOS (release:7.9.2009) with htcondor workflow management system.
If you fancy using MacOS or Windows, please design your set up accordingly.

## Clone github repository
```
cd ~/NatProt
git clone https://github.com/bartosovic-lab/nanoscope
```
> If you encounter authentication errors, you need to create a [personal access token](https://docs.github.com/fr/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)

## Prepare environment
A conda environment will be used to set up an isolated architecture reducing troubleshooting.
Conda can be installed via miniconda following [miniconda guidelines](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

### Create conda environment
```
conda create --name nanoscope_base -f nanoscope/envs/nanoscope_base.yaml 
```
Sit in your new environment till the end of the procedure.
```
conda activate nanoscope_base
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
fasterq-dump -f -e 1 --split-files --include-technical -o SRR18305888.fastq SRR18305888
fasterq-dump -f -e 1 --split-files --include-technical -o SRR18305889.fastq SRR18305889

mv SRR18305888_1.fastq sample_P23209_001_1001_S1_L001_I1_001.fastq
mv SRR18305888_2.fastq sample_P23209_001_1001_S1_L001_R1_001.fastq
mv SRR18305888_3.fastq sample_P23209_001_1001_S1_L001_R2_001.fastq
mv SRR18305888_4.fastq sample_P23209_001_1001_S1_L001_R3_001.fastq
mv SRR18305889_1.fastq sample_P23209_001_1001_S1_L002_I1_001.fastq
mv SRR18305889_2.fastq sample_P23209_001_1001_S1_L002_R1_001.fastq
mv SRR18305889_3.fastq sample_P23209_001_1001_S1_L002_R2_001.fastq
mv SRR18305889_4.fastq sample_P23209_001_1001_S1_L002_R3_001.fastq

gzip *.fastq
mkdir -p ./fastq/sample_P23209/
mv *.fastq.gz ./fastq/sample_P23209/
```

```
fasterq-dump -f -e 1 --split-files --include-technical -o SRR18305884.fastq SRR18305884
fasterq-dump -f -e 1 --split-files --include-technical -o SRR18305885.fastq SRR18305885

mv SRR18305884_1.fastq sample_P24004_002_1001_S1_L001_I1_001.fastq
mv SRR18305884_2.fastq sample_P24004_002_1001_S1_L001_R1_001.fastq
mv SRR18305884_3.fastq sample_P24004_002_1001_S1_L001_R2_001.fastq
mv SRR18305884_4.fastq sample_P24004_002_1001_S1_L001_R3_001.fastq
mv SRR18305885_1.fastq sample_P24004_002_1001_S1_L002_I1_001.fastq
mv SRR18305885_2.fastq sample_P24004_002_1001_S1_L002_R1_001.fastq
mv SRR18305885_3.fastq sample_P24004_002_1001_S1_L002_R2_001.fastq
mv SRR18305885_4.fastq sample_P24004_002_1001_S1_L002_R3_001.fastq

gzip *.fastq
mkdir -p ./fastq/sample_P24004/
mv *.fastq.gz ./fastq/sample_P24004/
```

> This step takes a while...



## Install cellranger
Cellranger is a bit more complex to fit inside a conda environment and is generally heavy to store.
On most HPCs running bioinformatic pipelines, Cellranger is already installed.
However if you wish to run the pipeline on a separated workstation, you can follow [10Xgenomics guidelines](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation).

> If cellranger is not installed on your favorite HPC, please contact your HPC support.

## HPC profiles
The pipeline is implemented in workflow management software known as snakemake.
It communicates with HPCs to run paralellized jobs to speed up the process.
Like previously mentionned at the beginning of the [Set Up](#set-up), the conda environment has been built to create a communication between snakemake and htcordor scheduler, therefore, htcondor package has been installed in the conda environment.

For htcondor workflow management, we will follow these [guidelines](https://github.com/Snakemake-Profiles/htcondor)
```
mkdir -p ~/.config/snakemake
cd ~/.config/snakemake

cookiecutter https://github.com/Snakemake-Profiles/htcondor.git
```

> At `profile_name [htcondor]` press `enter` and select a path for your log files, something like `~/condor_jobs`


For slurm workflow management, follow these [guidelines](https://github.com/Snakemake-Profiles/slurm)
```
mkdir -p ~/.config/snakemake
cd ~/.config/snakemake

template="gh:Snakemake-Profiles/slurm"
cookiecutter --output-dir "$profile_dir" "$template"
```


If your HPC is running on a different scheduler, you can install different style of profiles, like [slurm](https://github.com/Snakemake-Profiles/slurm).

> Do not forget to change the profile in your snakemake command line : 
```
--profile slurm
```
> If you do not have access to an HPC, you can simply remove the profile options from the snakemake command line and ignore this section.

## Changing parameters

### Creation of temporary folder
This temporary directory will be given to the config file to handle snakemake temporary outputs
```
mkdir -p ~/tmp/NatProt
```
### Modify config.yaml
The config.yaml file in the github repository will tailored the workflow according to its content.
Here is showned the config file for the [downloaded fastq files](#download-sra).

One can change the name of the samples as well as the path of the fastq files and the associated modalities. General information can also be tweaked, such as the temporary directory and conda environment used by snakemake, and parameters related to cellranger binary and reference location.

**NB**: this tutorial assumes the cellranger atac reference has already been generated. If not, please follow cellranger instructions [here](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/advanced/references)


Go back to the git directory
```
cd ~/NatProt/nanoscope
```

You can find the config.yaml file in `config/config.yaml`

```
samples:
  sample_P23209:
    fastq_path:
      ~/NatProt/Data/fastq/sample_P23209/
    barcodes:
      ATAC: TATAGCCT
      H3K27ac: ATAGAGGC
      H3K27me3: CCTATCCT

  sample_P24004:
    fastq_path:
      ~/NatProt/Data/fastq/sample_P24004/
    barcodes:
      ATAC: TATAGCCT
      H3K27ac: ATAGAGGC
      H3K27me3: CCTATCCT

general:
  tempdir: ~/tmp/NatProt
  cellranger_software: /data/bin/cellranger-atac
  cellranger_ref: /data/ref/cellranger-atac/refdata-cellranger-atac-mm10-2020-A-2.0.0/
```

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
cd ~/NatProt/nanoscope
snakemake --snakefile workflow/Snakefile_preprocess.smk --cores 16 --profile htcondor -p
```

# Downstream analysis
Once that the [Set Up](#set-up) and [Preprocessing](#preprocessing) steps are succcesfully completed, data is ready for downstream analysis. In this part of the documentation we provide a vignette on how to perform the downstream analyses, from raw data to identification and annotation of the different cell states in the dataset in analysis.

Please, follow the [vignette](https://fansalon.github.io/vignette_single-cell-nanoCT.html).
