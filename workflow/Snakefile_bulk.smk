include: 'Snakefile_prep2.smk'
include: 'Snakefile_nanoscope.smk'



rule all:
    input:
        bigwig       = expand('{sample}/{modality}/mapping_out/{modality}_merged.bw',sample = config['name'], modality = list(set(barcodes_dict.values()))),
        macs_3column = expand('{sample}/{modality}/peaks/macs2/{modality}_peaks.broadPeak',sample = config['name'], modality = list(set(barcodes_dict.values())))

rule trim_trim_galore:
    input:
        fastq_R1 = '{sample}/fastq_debarcoded/barcode_{barcode}/{prefix}_{number}_{lane}_R1_{suffix}',
        fastq_R2 = '{sample}/fastq_debarcoded/barcode_{barcode}/{prefix}_{number}_{lane}_R3_{suffix}',
    output:
        out_R1 = '{sample}/{modality}_{barcode}/fastq_trimmed/{prefix}_{number}_{lane}_R1_{suffix}',
        out_R2 = '{sample}/{modality}_{barcode}/fastq_trimmed/{prefix}_{number}_{lane}_R3_{suffix}',
    params:
        outdir = str(Path('{sample}/{modality}_{barcode}/fastq_trimmed/{prefix}_{number}_{lane}_R1_{suffix}').parents[0]),
        # OUT_R1 = '{sample}_{lane}_val_1'
        # OUT_R2 = '{prefix}_{number}_{lane}_R3_{suffix}',
    conda: "../envs/nanoscope_general.yaml"
    threads: 8
    resources:
        mem_mb=8000
    shell:
        """
        trim_galore --gzip --cores {threads} --fastqc --paired -o {params.outdir} --basename {wildcards.sample}_{wildcards.lane} {input.fastq_R1} {input.fastq_R2};
        TRIM_OUT_R1=`ls {params.outdir}/*_{wildcards.lane}*_val_1.fq.gz`;
        TRIM_OUT_R2=`ls {params.outdir}/*_{wildcards.lane}*_val_2.fq.gz`;
        mv $TRIM_OUT_R1 {output.out_R1}
        mv $TRIM_OUT_R2 {output.out_R2}
        """

rule map_bwa:
    input:
        read1 = '{sample}/{modality}_{barcode}/fastq_trimmed/{prefix}_{number}_{lane}_R1_{suffix}',
        read2 = '{sample}/{modality}_{barcode}/fastq_trimmed/{prefix}_{number}_{lane}_R3_{suffix}',
        index = str(Path(config['general']['cellranger_ref'] + '/fasta/genome.fa').resolve()),
    output:
        bam = '{sample}/{modality}_{barcode}/mapping_out/{prefix}_{number}_{lane}_{suffix}_mapped.bam',
    conda: "../envs/nanoscope_general.yaml"
    threads: 8
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
        '{sample}/{modality}_{barcode}/mapping_out/{prefix}_{number}_{lane}_{suffix}_mapped.bam',
    output:
        bam_sorted = '{sample}/{modality}_{barcode}/mapping_out/{prefix}_{number}_{lane}_{suffix}_sorted.bam',
        bam_index  = '{sample}/{modality}_{barcode}/mapping_out/{prefix}_{number}_{lane}_{suffix}_sorted.bam' + '.bai',
    threads: 8
    conda: "../envs/nanoscope_general.yaml"
    resources:
        mem_mb=32000
    shell:
        "samtools sort -o {output.bam_sorted} -@ {threads} {input} && samtools index {output.bam_sorted} "

rule merge_mapped:
  input:
    bam = lambda wildcards:merge_bam_inputs(config['fastq_path'],sample = wildcards.sample,modality=wildcards.modality,barcodes_dict = barcodes_dict)
  output:
    bam = '{sample}/{modality}/mapping_out/{modality}_merged.bam'
  threads: 8
  conda: "../envs/nanoscope_general.yaml"
  resources:
    mem_mb=32000
  shell:
    "samtools merge -@ {threads} {output.bam} {input.bam} && samtools index {output.bam}"

rule bam_to_bigwig:
  input:
    '{sample}/{modality}/mapping_out/{modality}_merged.bam'
  output:
    '{sample}/{modality}/mapping_out/{modality}_merged.bw'
  conda: "../envs/nanoscope_general.yaml"
  threads: 8
  resources:
    mem_mb=16000
  shell:
    "bamCoverage -b {input} -o {output} -p {threads} --normalizeUsing RPKM"

# rule run_macs_broad_bulk:
#     input:
#         bam_merged_wildcard
#     output:
#         macs_wildcard
#     params:
#         macs_outdir = str(Path(macs_wildcard).parents[0]),
#         macs_genome = config['general']['macs_genome']
#     conda: '../envs/nanoscope_general.yaml'
#     resources:
#         mem_mb = 16000
#     shell:
#         'macs2 callpeak -t {input} -g {params.macs_genome} -f BAMPE -n {wildcards.sample}_{wildcards.modality} '
#         '--outdir {params.macs_outdir} --llocal 100000 --keep-dup 1 --broad-cutoff 0.1 '
#         '--max-gap 1000 --broad 2>&1 '

rule macs_per_modality:
    input:
        '{sample}/{modality}/mapping_out/{modality}_merged.bam'
    output:
        '{sample}/{modality}/peaks/macs2/{modality}_peaks.broadPeak'
    params:
        macs_outdir = str(Path('{sample}/{modality}/peaks/macs2/{modality}_peaks.broadPeak').parents[0]),
        macs_genome = config['general']['macs_genome']
    conda: '../envs/nanoscope_general.yaml'
    resources:
        mem_mb = 32000
    shell:
        'macs2 callpeak -t {input} -g {params.macs_genome} -f BAMPE -n {wildcards.modality} '
        '--outdir {params.macs_outdir} --llocal 100000 --keep-dup 1 --broad-cutoff 0.1 '
        '--max-gap 1000 --broad 2>&1 '

rule create_fasta_index:
    output:
        fasta_index = '{sample}/general/{sample}_genome_index.fa.fai'
    params:
        fasta       = os.path.normpath(config['general']['cellranger_ref'] + '/fasta/genome.fa'),
    conda: '../envs/nanoscope_general.yaml'
    shell:
        'samtools faidx -o {output.fasta_index} {params.fasta}; '

rule peaks_to_3column_bed:
    input:
        peaks       = '{sample}/{modality}/peaks/macs2/{modality}_peaks.broadPeak',
        fasta_index = os.path.normpath(config['general']['cellranger_ref'] + '/fasta/genome.fa.fai'),
    output:
        '{sample}/{modality}/peaks/macs2/{modality}_peaks_3column.bed'
    conda: '../envs/nanoscope_general.yaml'
    shell:
        'cut -f1-3 {input.peaks} | bedtools sort -faidx {input.fasta_index} -i - > {output}'