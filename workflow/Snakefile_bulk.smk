include: 'Snakefile_prep.smk'
include: 'Snakefile_nanoscope.smk'

rule all:
    input:
        trim         = [run.samples[s].trimmed_fastq_all[i].path for s in run.samples_list for i,item in enumerate(run.samples[s].trimmed_fastq_all)],
        bowtie2      = [run.samples[s].bowtie2_bam_all for s in run.samples_list],
        bam_sorted   = [run.samples[s].bam_sorted_all for s in run.samples_list],
        bam_merged   = [run.samples[s].bam_merged_all for s in run.samples_list],
        bigwigs      = [run.samples[s].bigwig_all for s in run.samples_list],
        peaks        = [run.samples[s].macs_all for s in run.samples_list],
        peaks_merged = [run.samples[s].macs_merged_accross_modality_all for s in run.samples_list],

rule trim_trim_galore:
    input:
        fastq_R1 = lambda wildcards: run.samples[wildcards.sample].debarcoded_fastq_dict[wildcards.lane][wildcards.barcode]['R1'].path,
        fastq_R2 = lambda wildcards: run.samples[wildcards.sample].debarcoded_fastq_dict[wildcards.lane][wildcards.barcode]['R3'].path,
    output:
        out_R1 = temp(trimmed_fastq_output['R1']),
        out_R2 = temp(trimmed_fastq_output['R3']),
    params:
        outdir = str(Path(trimmed_fastq_wildcard).parents[0]),
        OUT_R1 = lambda wildcards: run.samples[wildcards.sample].trimmed_fastq_dict[wildcards.lane][wildcards.barcode]['R1'].path,
        OUT_R2 = lambda wildcards: run.samples[wildcards.sample].trimmed_fastq_dict[wildcards.lane][wildcards.barcode]['R3'].path
    conda: "../envs/bulk/nanoscope_trim.yaml"
    threads: 8
    resources:
        mem_mb=8000
    shell:
        """
        trim_galore --cores {threads} --fastqc --paired -o {params.outdir} --basename {wildcards.sample}_{wildcards.lane} {input.fastq_R1} {input.fastq_R2};
        TRIM_OUT_R1=`ls {params.outdir}/*_{wildcards.lane}*_val_1.fq.gz`;
        TRIM_OUT_R2=`ls {params.outdir}/*_{wildcards.lane}*_val_2.fq.gz`;
        mv $TRIM_OUT_R1 {params.OUT_R1}
        mv $TRIM_OUT_R2 {params.OUT_R2}
        """

rule map_bwa:
    input:
        read1 = lambda wildcards: run.samples[wildcards.sample].trimmed_fastq_dict[wildcards.lane][wildcards.barcode]['R1'].path,
        read2 = lambda wildcards: run.samples[wildcards.sample].trimmed_fastq_dict[wildcards.lane][wildcards.barcode]['R3'].path,
        index = bwa_index
    output:
        bam = temp(bowtie2_map_wildcard),
    conda: "../envs/bulk/nanoscope_bwa.yaml"
    threads: 20
    resources:
        mem_mb=32000
    shell:
        """
        bwa mem -t {threads} -M \
                {input.index} \
                {input.read1} \
                {input.read2} | samtools view -bSh > {output.bam}
        """


rule bam_sort_and_index:
    input:
        bowtie2_map_wildcard
    output:
        bam_sorted = temp(bam_sorted_wildcard),
        bam_index  = temp(bam_sorted_wildcard + '.bai')
    threads: 16
    conda: "../envs/nanoscope_samtools.yaml"
    resources:
        mem_mb=32000
    shell:
        "samtools sort -o {output.bam_sorted} -@ {threads} {input} && samtools index {output.bam_sorted} "

rule merge_mapped:
  input:
    bam = lambda wildcards: [bam_sorted_wildcard.replace('{lane}',l).replace('{barcode}',b) for l in run.samples[wildcards.sample].all_lanes for b in run.samples[wildcards.sample].reverse_barcodes_dict[wildcards.modality]],
  output:
    bam = bam_merged_wildcard
  threads: 16
  conda: "../envs/nanoscope_samtools.yaml"
  resources:
    mem_mb=32000
  shell:
    "samtools merge -@ {threads} {output.bam} {input.bam} && samtools index {output.bam}"

rule bam_to_bigwig:
  input:
    bam_merged_wildcard
  output:
    bigwig_wildcard
  conda: "../envs/nanoscope_deeptools.yaml"
  threads: 16
  resources:
    mem_mb=16000
  shell:
    "bamCoverage -b {input} -o {output} -p {threads} --normalizeUsing RPKM"

rule run_macs_broad_bulk:
    input:
        bam_merged_wildcard
    output:
        macs_wildcard
    params:
        macs_outdir = str(Path(macs_wildcard).parents[0]),
        macs_genome = config['general']['macs_genome']
    conda: '../envs/nanoscope_deeptools.yaml'
    resources:
        mem_mb = 16000
    shell:
        'macs2 callpeak -t {input} -g {params.macs_genome} -f BAMPE -n {wildcards.sample}_{wildcards.modality} '
        '--outdir {params.macs_outdir} --llocal 100000 --keep-dup 1 --broad-cutoff 0.1 '
        '--max-gap 1000 --broad 2>&1 '

rule macs_per_modality:
    input:
        lambda wildcards: [bam_merged_wildcard.format(sample = wildcards.sample, modality = run.samples[wildcards.sample].barcodes_dict[b], barcode = b) for b in run.samples[wildcards.sample].barcodes_list]
    output:
        macs_merged_per_modality_wildcard
    params:
        macs_outdir = str(Path(macs_merged_per_modality_wildcard).parents[0]),
        macs_genome = config['general']['macs_genome']
    conda: '../envs/nanoscope_deeptools.yaml'
    resources:
        mem_mb = 16000
    shell:
        'macs2 callpeak -t {input} -g {params.macs_genome} -f BAMPE -n {wildcards.modality} '
        '--outdir {params.macs_outdir} --llocal 100000 --keep-dup 1 --broad-cutoff 0.1 '
        '--max-gap 1000 --broad 2>&1 '

rule create_fasta_index:
    output:
        fasta_index = fasta_index_wildcard
    params:
        fasta       = config['general']['cellranger_ref'] + '/fasta/genome.fa',
    conda: '../envs/nanoscope_samtools.yaml'
    shell:
        'samtools faidx -o {output.fasta_index} {params.fasta}; '

rule peaks_to_3column_bed:
    input:
        peaks       = macs_merged_per_modality_wildcard,
        fasta_index = fasta_index_wildcard
    output:
        macs_merged_per_modality_wildcard + '_3column.bed'
    conda: '../envs/nanoscope_samtools.yaml'
    shell:
        'cut -f1-3 {input.peaks} | bedtools sort -faidx {input.fasta_index} -i - > {output}'