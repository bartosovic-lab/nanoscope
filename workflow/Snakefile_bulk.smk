include: 'Snakefile_prep.smk'

rule all:
    input:
        trim = [run.samples[s].trimmed_fastq_all.path for s in run.samples_list]


rule trim_trim_galore:
    input:
        fastq_R1 = run.samples[wildcards.sample].debarcoded_fastq_by_lane[wildcards.lane][wildcards.modality]['R1'],
        fastq_R2 = run.samples[wildcards.sample].debarcoded_fastq_by_lane[wildcards.lane][wildcards.modality]['R3'],
    output:
        out_R1 = trimmed_fastq_wildcard['R1'],
        out_R2 = trimmed_fastq_wildcard['R2']
    params:
        outdir="{sample}/{modality}_{barcode}/temp/trimming/"
    conda: "../envs/bulk/nanoscope_trim.yaml"
    threads: 8
    resources:
        mem_mb=8000
    shell:
        "trim_galore --cores {threads} --fastqc --paired -o {params.outdir} --basename {wildcards.sample}_{wildcards.lane} {input}"


#
# rule map_bowtie2:
#     input:
#         read1="out/{sample}/trimming/{sample}_{lane}_val_1.fq.gz",
#         read2="out/{sample}/trimming/{sample}_{lane}_val_2.fq.gz",
#     output:
#         bam=temp("out/{sample}/{sample}_{lane}_mapped.bam"),
#         log="out/logs/{sample}/{sample}_{lane}_bowtie2_map.log"
#     params:
#         index=config['general']['bowtie2_index']
#     conda: "../envs/bulkCT_map.yaml"
#     threads: 16
#     resources:
#         mem_mb=32000
#     shell:
#         """
#         bowtie2 --threads {threads} \
#                 --dovetail \
#                 -x {params.index} \
#                 -1 {input.read1} \
#                 -2 {input.read2} 2> {output.log} | samtools view -bS > {output.bam}
#         """
#
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
