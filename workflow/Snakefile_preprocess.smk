include: 'Snakefile_prep.smk'

import sys
args = sys.argv
config_path = args[args.index("--configfile") + 1]


rule all_preprocess:
    input:
        cellranger=[
            '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'.format(sample=config['sample_name'],modality = m, barcode = b) for m,b in config['barcodes'].items()],
        # bigwig_all=[
        #     '{sample}/{modality}_{barcode}/bigwig/all_reads.bw'.format(sample=sample,modality=modality,barcode=
        #     barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        # macs_broad=[
        #     '{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak'.format(sample=sample,modality=modality,barcode=
        #     barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
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
        # cellranger_cleanup = [
        #     '{sample}/{modality}_{barcode}/_clean_cellranger.out'.format(sample=sample,modality=modality,barcode=
        #     barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        # matrix_peaks = [
        #     '{sample}/{modality}_{barcode}/matrix/matrix_peaks/'.format(sample=sample,modality=modality,barcode=
        #     barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        # matrix_bins = [
        #     '{sample}/{modality}_{barcode}/matrix/matrix_bin_{binsize}/'.format(sample=sample,modality=modality,barcode=
        #     barcodes_dict[sample][modality],binsize = binsize) for sample in samples_list for modality in barcodes_dict[sample].keys() for binsize in bins], 
        # matrix_genes = [
        #     '{sample}/{modality}_{barcode}/matrix/matrix_genes/'.format(sample=sample,modality=modality,barcode=
        #     barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        
        

rule split_R2:
    input:
        fastq=lambda wildcards: glob.glob(config['fastq_path'] + '/**/*{lane}*R[123]*.fastq.gz'.format(lane=wildcards.lane),recursive=True)
    output:
        '{sample}/debarcoded_fastq_temp/ATAC/{sample}_{number}_{lane}_R1_{suffix}',
        '{sample}/debarcoded_fastq_temp/ATAC/{sample}_{number}_{lane}_R2_{suffix}',
        '{sample}/debarcoded_fastq_temp/ATAC/{sample}_{number}_{lane}_R3_{suffix}',
        '{sample}/debarcoded_fastq_temp/nonATAC/{sample}_{number}_{lane}_R1_{suffix}',
        '{sample}/debarcoded_fastq_temp/nonATAC/{sample}_{number}_{lane}_R2_{suffix}',
        '{sample}/debarcoded_fastq_temp/nonATAC/{sample}_{number}_{lane}_R3_{suffix}',
        '{sample}/debarcoded_fastq_temp/nonATAC/{sample}_{number}_{lane}_R2_{suffix}_modality.gz',
        '{sample}/debarcoded_fastq_temp/noMatch/{sample}_{number}_{lane}_R1_{suffix}',
        '{sample}/debarcoded_fastq_temp/noMatch/{sample}_{number}_{lane}_R2_{suffix}',
        '{sample}/debarcoded_fastq_temp/noMatch/{sample}_{number}_{lane}_R3_{suffix}',
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 480
    params:
        script = workflow.basedir + '/scripts/split_ATAC.py',
    conda: '../envs/nanoscope_general.yaml'
    shell:
        """
        mkdir -p {wildcards.sample}/debarcoded_fastq_temp/;
        cd {wildcards.sample}/debarcoded_fastq_temp/;
        python3 {params.script} -i {input.fastq} --ATAC ATAC --nonATAC nonATAC --noMatch noMatch
        """

rule move_ATAC:
    input:
        R1 = '{sample}/debarcoded_fastq_temp/ATAC/{sample}_{number}_{lane}_R1_{suffix}',
        R2 = '{sample}/debarcoded_fastq_temp/ATAC/{sample}_{number}_{lane}_R2_{suffix}',
        R3 = '{sample}/debarcoded_fastq_temp/ATAC/{sample}_{number}_{lane}_R3_{suffix}',
    output:
        R1 = '{sample}/debarcoded_fastq/barcode_MeA/{sample}_{number}_{lane}_R1_{suffix}',
        R2 = '{sample}/debarcoded_fastq/barcode_MeA/{sample}_{number}_{lane}_R2_{suffix}',
        R3 = '{sample}/debarcoded_fastq/barcode_MeA/{sample}_{number}_{lane}_R3_{suffix}',
    shell:
        "mv {input.R1} {output.R1}; mv {input.R2} {output.R2}; mv {input.R3} {output.R3}"

rule config2barcodes:
    input:
        script=workflow.basedir + '/scripts/fastq_multx_barcodes_from_config.py',
    output:
        '{sample}/debarcoded_fastq_temp/{sample}_{lane}_barcodes.txt',
    params:
        config = config_path,
    threads: 1
    resources:
        mem_mb = 2000,
    conda: '../envs/nanoscope_general.yaml'
    shell:
        """
        python3 {input.script} --config {params.config} --output {output} --reverse-complement
        """

rule debarcode_fastq_multx:
    input:
        barcodes      = '{sample}/debarcoded_fastq_temp/{sample}_{lane}_barcodes.txt',
        R1            = '{sample}/debarcoded_fastq_temp/nonATAC/{sample}_{number}_{lane}_R1_{suffix}',
        R2_singlecell = '{sample}/debarcoded_fastq_temp/nonATAC/{sample}_{number}_{lane}_R2_{suffix}',
        R3            = '{sample}/debarcoded_fastq_temp/nonATAC/{sample}_{number}_{lane}_R3_{suffix}',
        R2_modality   = '{sample}/debarcoded_fastq_temp/nonATAC/{sample}_{number}_{lane}_R2_{suffix}_modality.gz',
    output: 
        R1  = expand('{{sample}}/debarcoded_fastq_temp/{{lane}}/barcode_{barcode}/{{sample}}_{{number}}_{{lane}}_{read}_{{suffix}}',barcode = [b for b in config['barcodes'].values() if b != 'MeA'], read=['R1','R2','R3'] ),
        #
    threads: 1
    params: 
        all_barcodes = [b for b in config['barcodes'].values() if b != 'MeA'],
        log = '{sample}/debarcoded_fastq/{sample}_{number}_{lane}_{suffix}_debarcode.log',
    conda: "../envs/nanoscope_general.yaml"
    shell:
        """
        for b in {params.all_barcodes}; do
          mkdir -p {wildcards.sample}/debarcoded_fastq_temp/{wildcards.lane}/barcode_$b/;
          done;
        mkdir -p {wildcards.sample}/debarcoded_fastq_temp/{wildcards.lane}/barcode_unmatched/; 
        fastq-multx {input.barcodes} {input.R2_modality} {input.R1} {input.R2_singlecell} {input.R3} \
        -o {wildcards.sample}/debarcoded_fastq_temp/{wildcards.lane}/barcode_%/{wildcards.sample}_{wildcards.number}_{wildcards.lane}_modality_{wildcards.suffix} \
        {wildcards.sample}/debarcoded_fastq_temp/{wildcards.lane}/barcode_%/{wildcards.sample}_{wildcards.number}_{wildcards.lane}_R1_{wildcards.suffix} \
        {wildcards.sample}/debarcoded_fastq_temp/{wildcards.lane}/barcode_%/{wildcards.sample}_{wildcards.number}_{wildcards.lane}_R2_{wildcards.suffix} \
        {wildcards.sample}/debarcoded_fastq_temp/{wildcards.lane}/barcode_%/{wildcards.sample}_{wildcards.number}_{wildcards.lane}_R3_{wildcards.suffix} &> {params.log} ;
        """

rule move_nano_CT:
    input:
        R1            = '{sample}/debarcoded_fastq_temp/{lane}/barcode_{barcode}/{sample}_{number}_{lane}_R1_{suffix}',
        R2            = '{sample}/debarcoded_fastq_temp/{lane}/barcode_{barcode}/{sample}_{number}_{lane}_R2_{suffix}',
        R3            = '{sample}/debarcoded_fastq_temp/{lane}/barcode_{barcode}/{sample}_{number}_{lane}_R3_{suffix}',
    output:
        R1 = '{sample}/debarcoded_fastq/barcode_{barcode}/{sample}_{number}_{lane}_R1_{suffix}',
        R2 = '{sample}/debarcoded_fastq/barcode_{barcode}/{sample}_{number}_{lane}_R2_{suffix}',
        R3 = '{sample}/debarcoded_fastq/barcode_{barcode}/{sample}_{number}_{lane}_R3_{suffix}',
    threads: 1
    shell:
        "mv {input.R1} {output.R1}; mv {input.R2} {output.R2}; mv {input.R3} {output.R3}"
      

rule run_cellranger:
    input:
        lambda wildcards: get_fastq_for_cellranger(config['fastq_path'], wildcards.sample, wildcards.modality, wildcards.barcode)
    output:
        bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        frag='{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        meta='{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv',
        peaks='{sample}/{modality}_{barcode}/cellranger/outs/peaks.bed',
    params:
        cellranger_software=config['general']['cellranger_software'],
        cellranger_ref=config['general']['cellranger_ref'],
        fastq_folder=lambda wildcards: os.getcwd() + '/{sample}/debarcoded_fastq/barcode_{barcode}'.format(sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode)
    threads: 20
    resources:
        mem_mb = 32000
    shell:
        # Clean up the temp debarcoded fastq folder
        'rm -rf {wildcards.sample}/debarcoded_fastq_temp/; '
        # Clean up previous cellranger run
        'rm -rf {wildcards.sample}/{wildcards.modality}_{wildcards.barcode}/cellranger/; '
        # Cd into the folder
        'cd {wildcards.sample}/{wildcards.modality}_{wildcards.barcode}/; '
        # Rum cellranger 
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
    resources:
        mem_mb = 32000
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
        lambda wildcards: ['{sample}/{modality}_{barcode}/cell_picking/metadata.csv'.format(sample = wildcards.sample, modality = m, barcode = b) for m,b in config['barcodes'].items()],
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
    

rule create_genebody_and_promoter_matrix:
    input:
        cellranger_gtf = config['general']['cellranger_ref'] + 'genes/genes.gtf.gz',
        frag           = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        cells          = '{sample}/all_cells.txt',
    output:
        genebody_gtf = '{sample}/{modality}_{barcode}/matrix/matrix_genes/genebody_and_promoter.gtf',
        genebody_bed = '{sample}/{modality}_{barcode}/matrix/matrix_genes/genebody_and_promoter.bed',
        gene_names   = '{sample}/{modality}_{barcode}/matrix/matrix_genes/gene_names.txt',
        features     = '{sample}/{modality}_{barcode}/matrix/matrix_genes/features.tsv.gz',
        matrix       = '{sample}/{modality}_{barcode}/matrix/matrix_genes/matrix.mtx.gz',
        barcodes     = '{sample}/{modality}_{barcode}/matrix/matrix_genes/barcodes.tsv',
        folder       = directory('{sample}/{modality}_{barcode}/matrix/matrix_genes/'),
    conda: '../envs/nanoscope_general.yaml'
    params:
        script = workflow.basedir + '/scripts/filter_cellranger_gtf_file.py',
    shell:
        'python3 {params.script} -i {input.cellranger_gtf} -o {output.genebody_gtf} -n {output.gene_names};'
        'cut -f 1,4,5 {output.genebody_gtf} > {output.genebody_bed};'
        'fragtk matrix -f {input.frag} -b {output.genebody_bed} -c {input.cells} -o {output.folder}; '
