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
        peaks_overlap=[
            '{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        barcodes_sum=[
            '{sample}/{modality}_{barcode}/barcode_metrics/all_barcodes.txt'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        cell_pick=[
            '{sample}/{modality}_{barcode}/cell_picking/metadata.csv'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        noLA_bam=[
            '{sample}/{modality}_{barcode}/cellranger/outs/fragments_noLA_duplicates.tsv.gz'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        cellranger_cleanup = [
            '{sample}/{modality}_{barcode}/_clean_cellranger.out'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        all_cells = expand('{sample}/all_cells.txt', sample=samples_list),
        matrix = [generate_matrix_out(sample = sample, modality = modality, barcode = barcodes_dict[sample][modality]) 
            for sample in samples_list for modality in barcodes_dict[sample].keys()],


rule demultiplex:
    input:
        script=workflow.basedir + '/scripts/debarcode.py',
        fastq=lambda wildcards: glob.glob(config['samples'][wildcards.sample]['fastq_path'] + '/**/*{lane}*R[123]*.fastq.gz'.format(lane=wildcards.lane),recursive=True)
    output:
        '{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/{sample}_{number}_{lane}_R1_{suffix}',
        '{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/{sample}_{number}_{lane}_R2_{suffix}',
        '{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/{sample}_{number}_{lane}_R3_{suffix}',
    params:
        nbarcodes=lambda wildcards: len(config['samples'][wildcards.sample]['barcodes']),
        out_folder=lambda wildcards: '{sample}/{modality}_{barcode}/fastq/'.format(sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode),
    conda: '../envs/nanoscope_general.yaml'
    shell:
        "python3 {input.script} -i {input.fastq} -o {params.out_folder} --single_cell --barcode {wildcards.barcode} --name {wildcards.sample} 2>&1"

rule run_cellranger:
    input:
        lambda wildcards: get_fastq_for_cellranger(config['samples'][wildcards.sample]['fastq_path'],sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode)
    output:
        bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        frag='{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        meta='{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv',
        peaks='{sample}/{modality}_{barcode}/cellranger/outs/peaks.bed',
    params:
        cellranger_software=config['general']['cellranger_software'],
        cellranger_ref=config['general']['cellranger_ref'],
        fastq_folder=lambda wildcards: os.getcwd() + '/{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/'.format(sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode)
    threads: 20
    resources:
        mem_mb = 32000
    shell:
        'rm -rf {wildcards.sample}/{wildcards.modality}_{wildcards.barcode}/cellranger/; '
        'cd {wildcards.sample}/{wildcards.modality}_{wildcards.barcode}/; '
        '{params.cellranger_software} count --id cellranger --reference {params.cellranger_ref} --fastqs {params.fastq_folder}'

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
        bam   = temp('{sample}/{modality}_{barcode}/cellranger/outs/possorted_noLA_duplicates_bam.bam'),
        index = temp('{sample}/{modality}_{barcode}/cellranger/outs/possorted_noLA_duplicates_bam.bam.bai'),
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

rule bam_to_bw: # For QC reasons
    input:
        cellranger_bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'
    output:
        bigwig='{sample}/{modality}_{barcode}/bigwig/all_reads.bw'
    threads: 20
    conda: '../envs/nanoscope_general.yaml'
    resources:
        mem_mb = 16000
    shell:
        'bamCoverage -b {input.cellranger_bam} -o {output.bigwig} -p {threads} --minMappingQuality 5 '
        ' --binSize 50 --centerReads --smoothLength 250 --normalizeUsing RPKM --ignoreDuplicates --extendReads'

rule run_macs_broad:
    input:
        cellranger_bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'
    output:
        broad_peaks='{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak'
    params:
        macs_outdir='{sample}/{modality}_{barcode}/peaks/macs_broad/',
        macs_genome=config['general']['macs_genome']
    conda: '../envs/nanoscope_general.yaml'
    resources:
        mem_mb = 16000
    shell:
        'macs2 callpeak -t {input} -g {params.macs_genome} -f BAMPE -n {wildcards.modality} '
        '--outdir {params.macs_outdir} --llocal 100000 --keep-dup 1 --broad-cutoff 0.1 '
        '--max-gap 1000 --broad 2>&1 '

rule barcode_metrics_peaks:
    input:
        bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        peaks='{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak',
    output:
        overlap='{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt'
    params:
        get_cell_barcode=workflow.basedir + '/scripts/get_cell_barcode.awk',
        add_sample_to_list=workflow.basedir + '/scripts/add_sample_to_list.py',
        tmpdir=config['general']['tempdir']
    conda: '../envs/nanoscope_general.yaml'
    resources:
        mem_mb = 16000
    shell:
        'bedtools intersect -abam {input.bam} -b {input.peaks} -u | samtools view -f2 | '
        'awk -f {params.get_cell_barcode} | sed "s/CB:Z://g" |  '
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
        bcd_all='{sample}/{modality}_{barcode}/barcode_metrics/all_barcodes.txt',
        bcd_peak='{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt',
        peaks='{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak',
        metadata='{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv',
        fragments='{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
    output:
        '{sample}/{modality}_{barcode}/cell_picking/cells_10x.png',
        '{sample}/{modality}_{barcode}/cell_picking/cells_picked.png',
        '{sample}/{modality}_{barcode}/cell_picking/metadata.csv',
    params:
        script=workflow.basedir + '/scripts/pick_cells.R',
        out_prefix='{sample}/{modality}_{barcode}/cell_picking/',
    resources:
        mem_mb = 25000
    conda: '../envs/nanoscope_general.yaml'
    shell:
        "Rscript {params.script} --metadata {input.metadata} --fragments {input.fragments} --bcd_all {input.bcd_all} --bcd_peak {input.bcd_peak} --modality {wildcards.modality} --sample {wildcards.sample} --out_prefix {params.out_prefix}"

rule clean_cellranger_output:
    input:
        bam = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        frag = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        meta = '{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv',
        peaks = '{sample}/{modality}_{barcode}/cellranger/outs/peaks.bed',
    output:
        '{sample}/{modality}_{barcode}/_clean_cellranger.out'
    params:
        cellranger_folder = '{sample}/{modality}_{barcode}/cellranger/'
    shell:
        'touch {params.cellranger_folder}/tmp.txt;'                     # Create temp empty file to avoid error if the directory is empty
        'ls -d  {params.cellranger_folder}/* | grep -v outs | xargs rm -r; '
        'touch {output}'

rule get_cells:
    input:
        lambda wildcards: ['{sample}/{modality}_{barcode}/cell_picking/metadata.csv'.format(sample = wildcards.sample, 
                                                                                            modality = modality, 
                                                                                            barcode = barcodes_dict[wildcards.sample][modality]) 
            for modality in barcodes_dict[wildcards.sample].keys()]
    output:
        cells = '{sample}/all_cells.txt'
    params:
        script = workflow.basedir + '/scripts/get_passed_cells_barcodes.awk'
    shell:
        "cat {input} | awk -f {params.script} | sed 's/\"//g'| sort | uniq > {output.cells}"

rule create_matrix_peaks:
    # Requires fragtk installation (https://github.com/stuart-lab/fragtk)
    input:
        frag  = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        peaks = '{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak',
        cells = '{sample}/all_cells.txt'
    output:
        features = '{sample}/{modality}_{barcode}/matrix/matrix_peaks/features.tsv.gz',
        matrix   = '{sample}/{modality}_{barcode}/matrix/matrix_peaks/matrix.mtx.gz',
        barcodes = '{sample}/{modality}_{barcode}/matrix/matrix_peaks/barcodes.tsv',
        folder   = directory('{sample}/{modality}_{barcode}/matrix/matrix_peaks/'),
    shell:
        'fragtk matrix -f {input.frag} -b {input.peaks} -c {input.cells} -o {output.folder}'

rule create_matrix_bins:
    input:
        bam   = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        frag  = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        cells = '{sample}/all_cells.txt'
    output:
        features   = '{sample}/{modality}_{barcode}/matrix/matrix_bin_{bins}/features.tsv.gz',
        matrix     = '{sample}/{modality}_{barcode}/matrix/matrix_bin_{bins}/matrix.mtx.gz',
        barcodes   = '{sample}/{modality}_{barcode}/matrix/matrix_bin_{bins}/barcodes.tsv',
        folder     = directory('{sample}/{modality}_{barcode}/matrix/matrix_bin_{bins}/'),
        chromsizes = temp('{sample}/{modality}_{barcode}/chromsizes_{bins}.txt'),
        windows    = temp('{sample}/{modality}_{barcode}/windows_{bins}.txt'),
    conda: '../envs/nanoscope_general.yaml'
    shell:
        """
        samtools idxstats {input.bam} | cut -f1-2 | awk '$2 != 0' > {output.chromsizes}; 
        bedtools makewindows -g {output.chromsizes} -w {wildcards.bins} > {output.windows}; 
        fragtk matrix -f {input.frag} -b {output.windows} -c {input.cells} -o {output.folder}; 
        """


