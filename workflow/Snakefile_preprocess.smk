include: 'Snakefile_prep.smk'

rule all_preprocess:
    input:
        cellranger=[
            '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        bigwig_all=[
            '{sample}/{modality}_{barcode}/bigwig/all_reads.bw'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        macs_broad=[
            '{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        fragments=[
            '{sample}/{modality}_{barcode}/fragments/fragments.tsv.gz'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        peaks_overlap=[
            '{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        barcodes_sum=[
            '{sample}/{modality}_{barcode}/barcode_metrics/all_barcodes.txt'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        cell_pick=[
            '{sample}/{modality}_{barcode}/cell_picking/metadata.csv'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],

rule demultiplex:
    input:
        script=workflow.basedir + '/scripts/debarcode.py',
        fastq=lambda wildcards: glob.glob(config['samples'][wildcards.sample]['fastq_path'] + '/**/*{lane}*R[123]*.fastq.gz'.format(lane=wildcards.lane),recursive=True)
    output:
        '{sample}/fastq/{modality}_{barcode}/barcode_{barcode}/{id}_{number}_{lane}_R1_{suffix}',
        '{sample}/fastq/{modality}_{barcode}/barcode_{barcode}/{id}_{number}_{lane}_R2_{suffix}',
        '{sample}/fastq/{modality}_{barcode}/barcode_{barcode}/{id}_{number}_{lane}_R3_{suffix}',
    params:
        nbarcodes=lambda wildcards: len(config['samples'][wildcards.sample]['barcodes']),
        out_folder=lambda wildcards: '{sample}/fastq/{modality}_{barcode}/'.format(sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode),
    conda:
        '../envs/debarcode.yaml'
    shell:
        "python3 {input.script} -i {input.fastq} -o {params.out_folder} --single_cell --barcode {wildcards.barcode} 2>&1"

rule run_cellranger:
    input:
        lambda wildcards: get_fastq_for_cellranger(config['samples'][wildcards.sample][
            'fastq_path'],sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode)
    output:
        bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        frag='{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        meta='{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv',
        peaks='{sample}/{modality}_{barcode}/cellranger/outs/peaks.bed',
    params:
        cellranger_software=config['general']['cellranger_software'],
        cellranger_ref=config['general']['cellranger_ref'],
        fastq_folder=lambda wildcards: os.getcwd() + '/{sample}/fastq/{modality}_{barcode}/barcode_{barcode}/'.format(sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode)
    threads: 40
    shell:
        'rm -rf {wildcards.sample}/cellranger/{wildcards.sample}_{wildcards.modality}_{wildcards.barcode}/; '
        'cd {wildcards.sample}/cellranger/; '
        '{params.cellranger_software} count --id {wildcards.sample}_{wildcards.modality}_{wildcards.barcode} --reference {params.cellranger_ref} --fastqs {params.fastq_folder}'

rule bam_to_bw:
    input:
        cellranger_bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'
    output:
        bigwig='{sample}/{modality}_{barcode}/bigwig/all_reads.bw'
    threads: 16
    shell:
        'bamCoverage -b {input.cellranger_bam} -o {output.bigwig} -p {threads} --minMappingQuality 5 '
        ' --binSize 50 --centerReads --smoothLength 250 --normalizeUsing RPKM --ignoreDuplicates --extendReads'

rule run_macs_broad:
    input:
        cellranger_bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'
    output:
        broad_peaks='{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak'
    params:
        macs_outdir='{sample}/{modality}_{barcode}/peaks/macs_broad/'
    shell:
        'macs2 callpeak -t {input} -g mm -f BAMPE -n {wildcards.modality} '
        '--outdir {params.macs_outdir} --llocal 100000 --keep-dup=1 --broad-cutoff=0.1 '
        '--min-length 1000 --max-gap 1000 --broad 2>&1 '

rule add_barcode_fragments:
    input:
        fragments='{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz'
    output:
        fragments='{sample}/{modality}_{barcode}/fragments/fragments.tsv.gz',
        index='{sample}/{modality}_{barcode}/fragments/fragments.tsv.gz.tbi',
    params:
        script=workflow.basedir + '/scripts/add_sample_to_fragments.py',
    shell:
        'python3 {params.script} {input.fragments} {wildcards.sample} | bgzip > {output.fragments}; '
        'tabix -p bed {output.fragments}'

rule barcode_overlap_peaks:
    input:
        bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        peaks='{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak',
    output:
        overlap='{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt'
    params:
        get_cell_barcode=workflow.basedir + '/scripts/get_cell_barcode.awk',
        add_sample_to_list=workflow.basedir + '/scripts/add_sample_to_list.py',
        tmpdir=config['general']['tempdir']
    shell:
        'bedtools intersect -abam {input.bam} -b {input.peaks} -u | samtools view -f2 | '
        'awk -f {params.get_cell_barcode} | sed "s/CB:Z://g" | python3 {params.add_sample_to_list} {wildcards.sample} | '
        'sort -T {params.tmpdir} | uniq -c > {output.overlap} && [[ -s {output.overlap} ]] ; '

rule barcode_metrics_all:
    input:
        bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
    output:
        all_bcd='{sample}/{modality}_{barcode}/barcode_metrics/all_barcodes.txt'
    params:
        get_cell_barcode=workflow.basedir + '/scripts/get_cell_barcode.awk',
        add_sample_to_list=workflow.basedir + '/scripts/add_sample_to_list.py',
        tmpdir=config['general']['tempdir']
    shell:
        ' samtools view -f2 {input.bam}| '
        'awk -f {params.get_cell_barcode} | sed "s/CB:Z://g" | python3 {params.add_sample_to_list} {wildcards.sample} | '
        'sort -T {params.tmpdir} | uniq -c > {output.all_bcd} && [[ -s {output.all_bcd} ]] ; '

rule cell_selection:
    input:
        bcd_all='{sample}/{modality}_{barcode}/barcode_metrics/all_barcodes.txt',
        bcd_peak='{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt',
        peaks='{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak',
        metadata='{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv',
        fragments='{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
    output:
        '{sample}/{modality}_{barcode}/cell_picking/cells_10x.png',
        '{sample}/{modality}_{barcode}/cell_picking/cells_picked.png',
        '{sample}/{modality}_{barcode}/cell_picking/cells_picked.bw',
        '{sample}/{modality}_{barcode}/cell_picking/cells_not_picked.bw',
        '{sample}/{modality}_{barcode}/cell_picking/metadata.csv',
    params:
        script=workflow.basedir + '/scripts/pick_cells.R',
        out_prefix='{sample}/{modality}_{barcode}/cell_picking/',
    shell:
        "Rscript {params.script} --metadata {input.metadata} --fragments {input.fragments} --bcd_all {input.bcd_all} --bcd_peak {input.bcd_peak} --modality {wildcards.modality} --sample {wildcards.sample} --out_prefix {params.out_prefix}"
