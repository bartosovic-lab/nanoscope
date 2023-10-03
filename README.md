# Nanoscope 
# Analysis pipeline for single-cell nano-CUT&Tag data
### Federico Ansaloni, Bastien Hervé, Marek Bartosovic
<hr>

# Overview
Nanoscope is a snakemake-based pipeline for analysis of single-cell nano-CT data. For more details on sc-nano-CT please check out the following publications:

- https://www.nature.com/articles/s41587-022-01535-4
  
- TODO: ADD Nature Protocols link here
------------------------
The pipeline is composed of three majors axes, [Set Up](#set-up), [Preprocessing](#preprocessing) and [Downstream Analysis](#downstream-analysis) which will be broken down into sub-units to facilitate the go-through.

# Quick use (TL;DR)

If you wish to quickly use the pipeline follow these steps:

```
mkdir my_new_project
cd my_new_project

# Clone the git repo
git clone https://github.com/bartosovic-lab/nanoscope

# Create conda environment
conda config --set channel_priority flexible
conda env create -f nanoscope/envs/nanoscope_base.yaml 
conda activate nanoscope_base

# Change file paths and other parameters in pipeline config file 
vim nanoscope/config/config.yaml

# Run the pipeline
snakemake --snakefile nanoscope/workflow/Snakefile_preprocess.smk --configfile nanoscope/config/config.yaml --cores 20 --jobs 20 -p --use-conda --rerun-incomplete --profile slurm
```
# Hardware and software requirements

High Performance Computing cluster
-   20 cores
- 128GB RAM

The whole project has been tested on 
- HPC linux cluster under CentOS (release:7.9.2009) with htcondor workflow management system 
- HPC linux cluster under CentOS (release:7.9.2009) with slurm workflow management system


If you fancy using MacOS, Windows or another linux distro, please design your set up accordingly.

# Set up
Timing: 30 minutes

### Clone github repository
```
mkdir ~/nanoCT_project
cd ~/nanoCT_project
git clone https://github.com/bartosovic-lab/nanoscope
```

### Prepare environment
A conda environment will be used to set up an isolated architecture reducing troubleshooting.
Conda can be installed via miniconda following [miniconda guidelines](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

Createing complex conda environment from a yaml can tak >1h. For faster installation use mamba instead of conda.

### Create the environment
```
conda config --set channel_priority flexible
conda install -c conda-forge mamba
mamba env create -f nanoscope/envs/nanoscope_base.yaml 
```

Sit in your new environment till the end of the procedure.

```
conda activate nanoscope_base
```

### Install cellranger
Cellranger is a bit more complex to fit inside a conda environment and is generally heavy to store.
On most HPCs running bioinformatic pipelines, Cellranger is already installed.
However if you wish to run the pipeline on a separated workstation, you can follow [10Xgenomics guidelines](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation).

> If cellranger is not installed on your favorite HPC, please contact your HPC support.

Pre compiled references for cellranger-atac for some species are available [here](https://support.10xgenomics.com/single-cell-atac/software/downloads/latest)

To generate cellranger-atac reference please follow [here](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/advanced/references)

### HPC profiles
The pipeline is implemented in workflow management software known as snakemake.
It communicates with HPCs to run parallelized jobs to speed up the process.
Like previously mentioned at the beginning of the [Set Up](#set-up), the conda environment has been built to create a communication between snakemake and htcordor/slurm scheduler, therefore, htcondor package has been installed in the conda environment.

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
cookiecutter --output-dir ~/.config/snakemake "$template"
```

> If you do not have access to an HPC, you can simply remove the profile options from the snakemake command line and ignore this section.


# Data availability
All raw and processed files can be found as supplementary files in the [GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198467).

Processed `seurat.rds` object can also be downloaded from [GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198467) and used in the analysis from [Downstream Analysis](#downstream-analysis) or custom analysis.

> In this tutorial,  [GSM5949206](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198467) (SRR18305888 and SRR18305889) and [GSM5949208](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198467) (SRR18305884, SRR18305885) will be processed


# Download the raw data
Timing: 5h

Required storage space for raw data: 75GB

The raw data as fastq files can be downloaded throughout the [SRA-Toolkit](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump).

SRA tools can be installed using ```conda install -c bioconda sra-tools```

### Navigate to a new directory to store the data
```
mkdir -p ~/nanoCT_project/Data
cd ~/nanoCT_project/Data
```
### Download SRA
```
fasterq-dump -f -e 8 --split-files --include-technical -o SRR18305888.fastq SRR18305888
fasterq-dump -f -e 8 --split-files --include-technical -o SRR18305889.fastq SRR18305889

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
fasterq-dump -f -e 8 --split-files --include-technical -o SRR18305884.fastq SRR18305884
fasterq-dump -f -e 8 --split-files --include-technical -o SRR18305885.fastq SRR18305885

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

This step takes a while...

> Note: The file names should follow the standard illumina .fastq naming convention having for example the following format:

>SampleID_S1_L001_R1_001.fastq

> Otherwise the pipeline probably won't work. For more details on illumina fastq file naming please see: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm


# Pipeline parameters
Timing: 15 minutes

### Creation of temporary folder
This temporary directory will be given to the config file to handle snakemake temporary outputs
```
mkdir -p ~/tmp/nanoCT_project
```
### Modify config.yaml
The config.yaml file in the github repository will tailored the workflow according to its content.

Here is  the config file for the [downloaded fastq files](#download-sra).

One can change the name of the samples as well as the path of the fastq files and the associated modalities. General information can also be tweaked, such as the temporary directory and parameters related to cellranger binary and reference location.

**NB**: this tutorial assumes the cellranger-atac reference has already been generated. If not, please follow cellranger instructions [here](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/advanced/references) or download it directly from [cellranger website](https://support.10xgenomics.com/single-cell-atac/software/downloads/latest?).


Go back to the git directory
```
cd ~/nanoCT_project/nanoscope
```

You can find the config.yaml file in `config/config.yaml`

```
samples:
  sample_P23209:
    fastq_path:
      ~/nanoCT_project/Data/fastq/sample_P23209/
    barcodes:
      ATAC: TATAGCCT
      H3K27ac: ATAGAGGC
      H3K27me3: CCTATCCT

  sample_P24004:
    fastq_path:
      ~/nanoCT_project/Data/fastq/sample_P24004/
    barcodes:
      ATAC: TATAGCCT
      H3K27ac: ATAGAGGC
      H3K27me3: CCTATCCT

general:
  tempdir: ~/tmp/nanoCT_project
  cellranger_software: /data/bin/cellranger-atac
  cellranger_ref: /data/ref/cellranger-atac/refdata-cellranger-atac-mm10-2020-A-2.0.0/
  macs_genome: mm
```

# Preprocessing
Timing: ~8h on HPC - 16-20 cores per job and up to 100 jobs simultaneously running.

Required storage for processed data: ~150GB.

The first steps of processing from fastq files to cell picking will be done by the workflow management system, snakemake.
It will cover the following steps :
 - Modality Demultiplexing
 - Reads alignment with Cellranger
 - Create pseudobulk bigwig tracks for all cells (for QC purposes) 
 - Call peaks with Macs2
 - Create summary files with number of reads per cell (including PCR duplicates) and number of reads in peak regions
 - Redo cell selection, based on the summary files

All outputed files will be automatically generated and will be used to run the R markdown vignette below.

cd into project directory
```
cd ~/nanoCT_project/
```
Run snakemake with htcondor profile

```
snakemake --snakefile nanoscope/workflow/Snakefile_preprocess.smk --cores 16 --jobs 100 --profile htcondor -p --use-conda --configfile nanoscope/config/config.yaml
```
Slurm profile
```
snakemake --snakefile nanoscope/workflow/Snakefile_preprocess.smk --cores 16 --jobs 100 --profile slurm -p --use-conda --configfile nanoscope/config/config.yaml
```

In interactive shell (not recommended)
```
snakemake --snakefile nanoscope/workflow/Snakefile_preprocess.smk --cores 16 --jobs 100 -p --use-conda --configfile nanoscope/config/config.yaml
```


# Pipeline outputs

The important files generated by the pipeline are the following:

```
# Fragments file:
General:            $SAMPLE_ID/$MODALITY_$BARCODE/cellranger/outs/fragments.tsv.gz
Specific example:   sample_P23209/H3K27ac_ATAGAGGC/cellranger/outs/fragments.tsv.gz

# Fragments file without LA duplicates:
General:            $SAMPLE_ID/$MODALITY_$BARCODE/cellranger/outs/fragments_noLA_duplicates.tsv.gz
Specific example:   sample_P23209/H3K27ac_ATAGAGGC/cellranger/outs/fragments_noLA_duplicates.tsv.gz

# Peaks file:
General:            $SAMPLE_ID/$MODALITY_$BARCODE/peaks/macs_broad/$MODALITY_peaks.broadPeak
Specific example:   sample_P23209/H3K27ac_ATAGAGGC/peaks/macs_broad/H3K27ac_peaks.broadPeak

# Cell metadata
General:            $SAMPLE_ID/$MODALITY_$BARCODE/cell_picking/metadata.csv
Specific example:   sample_P23209/H3K27ac_ATAGAGGC/cell_picking/metadata.csv

# Cell picking results - column "passedMB" indicates cells that passed the cell picking algorithm
General:            $SAMPLE_ID/$MODALITY_$BARCODE/cell_picking/cells_picked.png
Specific example:   sample_P23209/H3K27ac_ATAGAGGC/cell_picking/cells_picked.png

# Bam file with mapped reads (possorted_bam.bam file from cellranger)
General:            $SAMPLE_ID/$MODALITY_$BARCODE/cellranger/outs/possorted_bam.bam
Specific example:   sample_P23209/H3K27ac_ATAGAGGC/cellranger/outs/possorted_bam.bam
```

- Cell  picking might be further optimised and other algorithms might be used/considered by the user. Cell picking will be under further development and changes might be introduced in the future.


# Downstream analysis
Once that the [Set Up](#set-up) and [Preprocessing](#preprocessing) steps are succcesfully completed, data is ready for downstream analysis. In this part of the documentation we provide a vignette on how to perform the downstream analyses, from raw data to identification and annotation of the different cell states in the dataset in analysis.

If you prefer to perform the analysis by using **peaks**, please, follow this [vignette](https://fansalon.github.io/vignette_single-cell-nanoCT.html).\
If you prefer to perform the analysis by using **bins**, please, follow this [vignette](https://fansalon.github.io/vignette_single-cell-nanoCT_bins.html).\
[Here](https://fansalon.github.io/comparison_vignette_single-cell-nanoCT.html) you can also find an additional vignette where the analyses run by using bins and peaks on the same dataset are compared.

The use of peaks or bins mostly depends on the workflow you are applying to your analysis. If in your analyses you need to integrate different datasets (ATAC-seq with nanoCT, or different nanoCT datasets), or different samples, then we suggest to use bins (bins are the same everywhere). A downside of using bins is that the signal is probably a bit less specific than when using peaks, and overall the clustering resolution could be a bit worse when using bins than peaks.


Note how the dataset analysed in this vignette is composed by 2 biological replicates of the same sample. Often experiments are designed in order to have multiple samples and/or biological conditions requiring data integration, differential analysis, etc. These are not part of this vignette, but will be implemented in other vignettes in the future.




## Difference between fragments.tsv vs fragments_noLA_duplicates.tsv comparison 
```
# Typically the level of LA duplicates is very low and does not need to be considered in the analysis

gunzip -c fragments.tsv.gz | wc -l
66284027

gunzip -c fragments_noLA_duplicates.tsv.gz | wc -l
64010202

# LA Duplication rate is (66,284,027 - 64,010,202)/66,284,027 = ~3.4 %

If you see highrer fraction of LA duplicates, consider using fragments_noLA_duplicates.tsv.gz file instead of fragments.tsv.gz in downstream analysis. 
```

