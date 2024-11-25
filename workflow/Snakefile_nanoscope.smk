include: 'Snakefile_prep.smk'
include: 'Snakefile_bulk.smk'

rule all_preprocess:
    input:
    # Single-cell outputs
        debarcode_out     = [run.samples[s].debarcoded_fastq_all[i].path for s in run.samples_list for i,item in enumerate(run.samples[s].debarcoded_fastq_all)],
        cellranger_out    = [run.samples[s].cellranger_all for s in run.samples_list],
        cell_picking_out  = [run.samples[s].cell_picking_all for s in run.samples_list],

    # Bulk outputs
        # bigwig_bulk = [run.samples[s].bigwig_all for s in run.samples_list],
        # peaks_overlap=[
        #     '{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt'.format(sample=sample,modality=modality,barcode=
        #     barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        # barcodes_sum=[
        #     '{sample}/{modality}_{barcode}/barcode_metrics/all_barcodes.txt'.format(sample=sample,modality=modality,barcode=
        #     barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        # cell_pick=[
        #     '{sample}/{modality}_{barcode}/cell_picking/metadata.csv'.format(sample=sample,modality=modality,barcode=
        #     barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        # noLA_bam=[
        #     '{sample}/{modality}_{barcode}/cellranger/outs/fragments_noLA_duplicates.tsv.gz'.format(sample=sample,modality=modality,barcode=
        #     barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],


rule demultiplex:
    input:
        script = workflow.basedir + '/scripts/debarcode.py',
        fastq_all = lambda wildcards: [run.samples[wildcards.sample].fastq_by_lane[wildcards.lane][r].path for r in run.samples[wildcards.sample].all_reads]
    output:
        debarcoded_fastq_all = [debarcoded_fastq_wildcard.format(sample = '{sample}', barcode = b, prefix = '{prefix}', number = '{number}', 
        lane = '{lane}', read = r, suffix = '{suffix}') for s in run.samples_list for b in run.samples[s].barcodes_dict.keys() for r in ['R1','R2','R3']]
    params:
        out_folder = str(Path(debarcoded_fastq_wildcard).parents[1]),
        all_barcodes = lambda wildcards: ' '.join(run.samples[wildcards.sample].barcodes_dict.keys()),
    conda: '../envs/nanoscope_general.yaml'
    threads: 1
    shell:
        "python {input.script} -i {input.fastq_all} -o {params.out_folder} --single_cell --barcode {params.all_barcodes}  2>&1"

rule run_cellranger:
    input:
        read1 = lambda wildcards: [run.samples[wildcards.sample].debarcoded_fastq_dict[l][wildcards.barcode]['R1'].path for l in run.samples[wildcards.sample].all_lanes],
        read2 = lambda wildcards: [run.samples[wildcards.sample].debarcoded_fastq_dict[l][wildcards.barcode]['R2'].path for l in run.samples[wildcards.sample].all_lanes],
        read3 = lambda wildcards: [run.samples[wildcards.sample].debarcoded_fastq_dict[l][wildcards.barcode]['R3'].path for l in run.samples[wildcards.sample].all_lanes],
        bed   = macs_merged_per_modality_wildcard + '_3column.bed'   # Peaks file comes from the Snakefile_bulk.smk pipeline
    output:
        fragments = cellranger_fragments_wildcard,
        bam       = cellranger_bam_wildcard,
        metadata  = cellranger_metadata_wildcard
    params:
        cellranger_software = config['general']['cellranger_software'],
        cellranger_ref      = config['general']['cellranger_ref'],
        fastq_folder        = str(Path(debarcoded_fastq_wildcard).parents[0].resolve()),
        bed_with_abspath    = str(Path(macs_merged_per_modality_wildcard + '_3column.bed').resolve())
    threads: 20
    resources:
        mem_mb = 32000
    shell:
        'rm -rf {wildcards.sample}/{wildcards.modality}_{wildcards.barcode}/cellranger/; '
        'cd {wildcards.sample}/{wildcards.modality}_{wildcards.barcode}/; '
        '{params.cellranger_software} count --id cellranger --reference {params.cellranger_ref} --fastqs {params.fastq_folder} --peaks={params.bed_with_abspath} --force-cells=5000'

rule barcode_metrics_peaks:
    input:
        bam   = cellranger_bam_wildcard,
        peaks = macs_merged_per_modality_wildcard + '_3column.bed'
    output:
        overlap_file = overlap_file_wildcard
    params:
        get_cell_barcode   = workflow.basedir + '/scripts/get_cell_barcode.awk',
        add_sample_to_list = workflow.basedir + '/scripts/add_sample_to_list.py',
        tmpdir             = config['general']['tempdir']
    conda: '../envs/nanoscope_general.yaml'
    resources:
        mem_mb = 16000
    shell:
        'bedtools intersect -abam {input.bam} -b {input.peaks} -u | samtools view -f2 | '
        'awk -f {params.get_cell_barcode} | sed "s/CB:Z://g" |  '
        'sort -T {params.tmpdir} | uniq -c > {output.overlap_file} && [[ -s {output.overlap_file} ]] ; '

rule barcode_metrics_all:
    input:
        bam  = cellranger_bam_wildcard
    output:
        all_bcd = bcd_stats_wildcard
    params:
        get_cell_barcode   = workflow.basedir + '/scripts/get_cell_barcode.awk',
        add_sample_to_list = workflow.basedir + '/scripts/add_sample_to_list.py',
        tmpdir             = config['general']['tempdir']
    conda: '../envs/nanoscope_general.yaml'
    resources:
        mem_mb = 16000
    shell:
        'mkdir -p {params.tmpdir}; '
        ' samtools view -f2 {input.bam}| '
        'awk -f {params.get_cell_barcode} | sed "s/CB:Z://g" |  '
        'sort -T {params.tmpdir} | uniq -c > {output.all_bcd} && [[ -s {output.all_bcd} ]] ; '



rule cell_selection:
    input:
        bcd_all   = bcd_stats_wildcard,
        bcd_peak  = overlap_file_wildcard,
        peaks     = macs_merged_per_modality_wildcard + '_3column.bed',
        metadata  = cellranger_metadata_wildcard,
        fragments = cellranger_fragments_wildcard
    output:
        cell_picking_cells_10x_wildcard,
        cell_picking_cells_nanoscope_wildcard,
        cell_picking_metadata_wildcard,
        # '{sample}/{modality}_{barcode}/cell_picking/cells_picked.bw',
        # '{sample}/{modality}_{barcode}/cell_picking/cells_not_picked.bw',
    params:
        script     = workflow.basedir + '/scripts/pick_cells.R',
        out_prefix = str(Path(cell_picking_metadata_wildcard).parents[0].resolve()), # This is the output folder
    resources:
        mem_mb = 25000
    conda: '../envs/nanoscope_general.yaml'
    shell:
        "Rscript {params.script} --metadata {input.metadata} --fragments {input.fragments} --bcd_all {input.bcd_all} --bcd_peak {input.bcd_peak} --modality {wildcards.modality} --sample {wildcards.sample} --out_prefix {params.out_prefix}"

############################################

rule cellranger_bam_to_namesorted:
    input:
        bam = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
    output:
        bam = temp('{sample}/{modality}_{barcode}/cellranger/outs/namesorted_bam.bam'),
    params:
        tempfolder = config['general']['tempdir']
    threads: 20
    resources:
        mem_mb = 32000
    conda: '../envs/nanoscope_general.yaml'
    shell:
        'samtools sort -T {params.tempfolder} -@ {threads} -n -o {output.bam} {input.bam} '

rule remove_LA_duplicates:
    input:
        bam = '{sample}/{modality}_{barcode}/cellranger/outs/namesorted_bam.bam',
    output:
        bam = temp('{sample}/{modality}_{barcode}/cellranger/outs/namesorted_noLA_duplicates_bam.bam'),
    params:
        script = workflow.basedir + '/scripts/remove_LA_duplicates.py',
    conda: '../envs/nanoscope_general.yaml'
    shell:
        'python3 {params.script} {input.bam} {output.bam}'

rule possort_noLA_bam_file:
    input:
        bam = '{sample}/{modality}_{barcode}/cellranger/outs/namesorted_noLA_duplicates_bam.bam',
    output:
        bam   = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_noLA_duplicates_bam.bam',
        index = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_noLA_duplicates_bam.bam.bai',
    conda: '../envs/nanoscope_general.yaml'
    threads: 20
    resources:
        mem_mb = 32000
    params:
        tempfolder = config['general']['tempdir']
    shell:
        'samtools sort -@ {threads} -T {params.tempfolder} -o {output.bam} {input.bam}; '
        'samtools index {output.bam} '

rule bam_noLA_to_fragments_noLA:
    input:
        bam = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_noLA_duplicates_bam.bam',
        index= '{sample}/{modality}_{barcode}/cellranger/outs/possorted_noLA_duplicates_bam.bam.bai',
    output:
        fragmments = temp('{sample}/{modality}_{barcode}/cellranger/outs/fragments_noLA_duplicates.tsv'),
    conda: '../envs/nanoscope_general.yaml'
    threads: 20
    resources:
        mem_mb=16000
    shell:
        'sinto fragments -b {input.bam} -f {output.fragmments} -p {threads}'

rule sort_sinto_output:
    input:
        fragments = '{sample}/{modality}_{barcode}/cellranger/outs/fragments_noLA_duplicates.tsv',
    output:
        fragments = '{sample}/{modality}_{barcode}/cellranger/outs/fragments_noLA_duplicates.tsv.gz',
        index      = '{sample}/{modality}_{barcode}/cellranger/outs/fragments_noLA_duplicates.tsv.gz.tbi',
    conda: '../envs/nanoscope_general.yaml'
    resources:
        mem_mb = 16000
    shell:
        'sort -k1,1 -k2,2n {input.fragments} | bgzip > {output.fragments}; '
        'tabix -p bed {output.fragments} '

# rule bam_to_bw: # For QC reasons
#     input:
#         cellranger_bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'
#     output:
#         bigwig='{sample}/{modality}_{barcode}/bigwig/all_reads.bw'
#     threads: 20
#     conda: '../envs/nanoscope_deeptools.yaml'
#     resources:
#         mem_mb = 16000
#     shell:
#         'bamCoverage -b {input.cellranger_bam} -o {output.bigwig} -p {threads} --minMappingQuality 5 '
#         ' --binSize 50 --centerReads --smoothLength 250 --normalizeUsing RPKM --ignoreDuplicates --extendReads'

# rule run_macs_broad:
#     input:
#         cellranger_bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'
#     output:
#         broad_peaks='{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak'
#     params:
#         macs_outdir='{sample}/{modality}_{barcode}/peaks/macs_broad/',
#         macs_genome=config['general']['macs_genome']
#     conda: '../envs/nanoscope_deeptools.yaml'
#     resources:
#         mem_mb = 16000
#     shell:
#         'macs2 callpeak -t {input} -g {params.macs_genome} -f BAMPE -n {wildcards.modality} '
#         '--outdir {params.macs_outdir} --llocal 100000 --keep-dup 1 --broad-cutoff 0.1 '
#         '--max-gap 1000 --broad 2>&1 '

