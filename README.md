# Nanoscope-mx
# Pipeline for analysis of multiplexed nano-CT data with overloading
### Marek Bartosovic
<hr>


Snakemake workflow for processing single-cell nano-CT libraries with 10x-style barcodes into:

1) demultiplexed FASTQs per barcode/modality

2) bulk alignments, and broad MACS2 peaks - this is later used for cellranger and cell calling

3) single-cell Cell Ranger ATAC outputs (BAM, fragments, singlecell.csv) using bulk peaks

4) QC and cell picking

5) count matrices (fixed genomic bins and gene body+promoter features)

## Quick start 

```
# 1) Clone
git clone https://github.com/bartosovic-lab/nanoscope-mx

# 2) (Recommended) use conda/mamba and let Snakemake manage envs

conda env create -f envs/nanoscope_base.yaml
conda activate nanoscope_base

# 3) Edit the config with your details

# 4) Run the pipeline 
mamba create -f env
snakemake -s Snakefile_nanoscope.smk \
  --use-conda --cores 16 \
  --configfile config/config.yaml
  ```

## Inputs 

FASTQ files

•	Libraries have three reads: R1, R2, R3 with read lengths 36, 48 and 36 respectively.

•	R2 contains barcodes - modality and single-cell barcode with the following arrangement.

It is assumed that sequencing is done using the nextera primer and bases 1-8 correspond to modality barcode and bases 30-45 correspond to the single-cell index. 

The barcodes are separated by linker sequence of fixed length TCGTCGGCAGCGTCTCCACGC

AATGATACGGCGACCACCGAGATCTACAC-NNNNNNNNNNNNNNNN-TCGTCGGCAGCGTCTCCACGC-NNNNNNNN-GCGATCGAGGACGGCAGATGTGTATAAGAGACAG

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;P5 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  sc-barcode  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Linker sequence  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  | &nbsp;&nbsp;Modality &nbsp;&nbsp;&nbsp;| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;       Mosaic end 



•	File naming assumes the standard illumina pattern:   
  
  \<PREFIX>_S\<NUMBER>_L\<LANE>_\<READ>_\<SUFFIX>.
  fastq.gz

•	If a library was re-sequenced, all files present in fastq folder will be used in analysis and reads merged for each modality barcode. e.g. 

``` 
# Single sequencing lane
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L001_R1_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L001_R2_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L001_R3_001.fastq.gz
```

``` 
# Two lanes - usual for NovaSeq X 1.5B
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L001_R1_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L001_R2_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L001_R3_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L002_R1_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L002_R2_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L002_R3_001.fastq.gz
```

``` 
# Multiple lanes from different sequencing runs
lrwxrwxrwx 1 mareba mareba 137 Oct 21 11:15 P33456_1001_S1_L001_I1_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 137 Oct 21 11:15 P33456_1001_S1_L001_R1_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 137 Oct 21 11:15 P33456_1001_S1_L001_R2_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 137 Oct 21 11:15 P33456_1001_S1_L001_R3_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 137 Oct 21 11:15 P33456_1001_S1_L002_I1_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 137 Oct 21 11:15 P33456_1001_S1_L002_R1_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 137 Oct 21 11:15 P33456_1001_S1_L002_R2_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 137 Oct 21 11:15 P33456_1001_S1_L002_R3_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L006_I1_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L006_R1_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L006_R2_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L006_R3_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L007_I1_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L007_R1_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L007_R2_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L007_R3_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L008_I1_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L008_R1_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L008_R2_001.fastq.gz
lrwxrwxrwx 1 mareba mareba 117 Oct 21 11:13 P37356_1014_S11_L008_R3_001.fastq.gz
```

In all cases - all the fastq files will be used in by the pipeline nad merged. If you wish to run independent samples, then put .fastq files in different folders and use multiple configs. 

## Pipeline config file 

Example minimal config

```{yaml}
# Required
name: "NanoCT_run1"               # sample/run name (used as top-level folder)
fastq_path: "/path/to/fastqs"     # root folder with FASTQs (recursive search)

barcodes:                         # barcode sequence => modality label
  ACGTAC: H3K27ac
  TGCATG: H3K4me3
  GATCTA: Input

general:
  cellranger_software: "cellranger-atac"    # path or name in $PATH
  cellranger_ref: "/refs/cellranger-atac/refdata"  # has fasta/genome.fa etc.
  macs_genome: "hs"                          # MACS2 -g (e.g., hs, mm, 2.7e9)
  tempdir: "/tmp"                            # for sorting/uniqs
  gtf: "https://.../gencode.vXX.annotation.gtf.gz" # used for gene bodies/promoters
  bin_sizes: [5000, 10000]                   # fixed bin matrices to generate
  ```

- name - arbitrary filesystem-safe 
- fastq_path - path to folder that contain .fastq.gz files from Illumina sequencing
- barcodes - mapping of barcodes to their respective modalities
- general - general parameters, like cellranger bin $PATH, cellranger ref path, macs genome (hs,mm,ce, dm) - see macs help, tempdir - path to temporary directory used in read sorting, gtf - url to gtf file containing gene annotations (gencode supported), bin_sizes - sizes of bins for which to construct the bin x cell matrix. 


## Outputs
```
NanoCT_run1/
  general/
    NanoCT_run1_barcodes_table.txt
    NanoCT_run1_genome_index.fa.fai
  fastq_split_R2/                          # R2 split into modality/singlecell
  fastq_debarcoded/
    barcode_<BARCODE>/                     # demultiplexed FASTQs per barcode
      <PREFIX>_S##_L###_R1_...
      <PREFIX>_S##_L###_R2_...             # single-cell R2 part
      <PREFIX>_S##_L###_R3_...
  <MODALITY>/                              # bulk per modality
    mapping_out/
      <MODALITY>_merged.bam
      <MODALITY>_merged.bam.bai
      <MODALITY>_merged.bw
    peaks/macs2/
      <MODALITY>_peaks.broadPeak
      <MODALITY>_peaks_3column.bed
  <MODALITY>_<BARCODE>/                    # single-cell per (modality, barcode)
    cellranger/outs/
      possorted_bam.bam
      fragments.tsv.gz
      singlecell.csv
      fragments_noLA_duplicates.tsv.gz(.tbi)   # after LA duplicate removal
      possorted_noLA_duplicates_bam.bam(.bai)
    barcode_metrics/
      peaks_barcodes.txt                    # reads per cell overlapping MACS peaks
      all_barcodes.txt                      # reads per cell overall
    cell_picking/
      cells_10x.png                         # Cell picking scatterplot using cellranger default
      cells_nanoscope.png                   # Cell picking scatterplot using nanoscope picker  
      metadata.csv                          # selected cells - replaces 10x cellranger metadata
    matrix/
      all_cells.txt
      matrix_bin_5000/ (or configured sizes)
        features.tsv.gz  barcodes.tsv  matrix.mtx.gz
      matrix_genebody_promoter/
        features.tsv.gz  barcodes.tsv  matrix.mtx.gz
```
## Pipeline scheme
!["Flowchart diagram showing the Nanoscope-mx pipeline workflow. Starting with multiplexed FASTQ files at the top, the process flows downward through three main stages. First stage: demultiplexing of reads by modality barcodes. Second stage: parallel processing of bulk alignments and peak calling for each modality. Third stage: single-cell processing including Cell Ranger ATAC analysis, cell calling, and generation of count matrices. Arrows connect the stages indicating data flow. Technical steps include BAM file generation, MACS2 peak calling, and creation of bin and gene body matrices."](img/image.png)
