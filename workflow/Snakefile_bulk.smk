include: 'Snakefile_prep.smk'
include: 'Snakefile_nanoscope.smk'

rule all:
    input:
        trim    = [run.samples[s].trimmed_fastq_all[i].path for s in run.samples_list for i,item in enumerate(run.samples[s].trimmed_fastq_all)],
        bowtie2 = [run.samples[s].bowtie2_bam_all for s in run.samples_list]


rule trim_trim_galore:
    input:
        fastq_R1 = lambda wildcards: run.samples[wildcards.sample].debarcoded_fastq_dict[wildcards.lane][wildcards.modality]['R1'].path,
        fastq_R2 = lambda wildcards: run.samples[wildcards.sample].debarcoded_fastq_dict[wildcards.lane][wildcards.modality]['R3'].path,
    output:
        out_R1 = trimmed_fastq_output['R1'],
        out_R2 = trimmed_fastq_output['R3'],
    params:
        outdir = trim_params_outdir,
        OUT_R1 = lambda wildcards: run.samples[wildcards.sample].trimmed_fastq_dict[wildcards.lane][wildcards.modality]['R1'].path,
        OUT_R2 = lambda wildcards: run.samples[wildcards.sample].trimmed_fastq_dict[wildcards.lane][wildcards.modality]['R3'].path
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

rule map_bowtie2:
    input:
        read1 = lambda wildcards: run.samples[wildcards.sample].trimmed_fastq_dict[wildcards.lane][wildcards.modality]['R1'].path,
        read2 = lambda wildcards: run.samples[wildcards.sample].trimmed_fastq_dict[wildcards.lane][wildcards.modality]['R3'].path
    output:
        bam = bowtie2_map_wildcard,
        log = bowtie2_map_wildcard.replace('.bam','.log')
    params:
        index = config['general']['bowtie2_index']
    conda: "../envs/bulk/nanoscope_map.yaml"
    threads: 16
    resources:
        mem_mb=32000
    shell:
        """
        bowtie2 --threads {threads} \
                --dovetail \
                -x {params.index} \
                -1 {input.read1} \
                -2 {input.read2} 2> {output.log} | samtools view -bS > {output.bam}
        """

# rule bam_sort_and_index:
#     input:
#         "out/{sample}/{sample}_{lane}_mapped.bam",
#     output:
#         bam_sorted=temp("out/{sample}/{sample}_{lane}_sorted.bam"),
#         bam_index=temp("out/{sample}/{sample}_{lane}_sorted.bam.bai"),
#     threads: 16
#     conda: "../envs/bulkCT_map.yaml"
#     resources:
#         mem_mb=32000
#     shell:
#         "samtools sort -o {output.bam_sorted} -@ {threads} {input} && samtools index {output.bam_sorted} "
