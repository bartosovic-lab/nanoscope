include: 'Snakefile_prep2.smk'
include: 'Snakefile_bulk.smk'

# barcodes_table_wildcard = '{sample}/general/{sample}_barcodes_table.txt'
# cellranger_fragments_wildcard         = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz'
# cellranger_bam_wildcard               = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'
# cellranger_metadata_wildcard          = '{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv'
# cellranger_cleanup_wildcard           = '{sample}/{modality}_{barcode}/cellranger/outs/_cellranger_cleanup'

# Single-cell wildcards
cellranger_fragments_wildcard         = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz'
cellranger_bam_wildcard               = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'
cellranger_metadata_wildcard          = '{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv'
cellranger_cleanup_wildcard           = '{sample}/{modality}_{barcode}/cellranger/outs/_cellranger_cleanup'
overlap_file_wildcard                 = '{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt'
bcd_stats_wildcard                    = '{sample}/{modality}_{barcode}/barcode_metrics/all_barcodes.txt'
cell_picking_cells_10x_wildcard       = '{sample}/{modality}_{barcode}/cell_picking/cells_10x.png'
cell_picking_cells_nanoscope_wildcard = '{sample}/{modality}_{barcode}/cell_picking/cells_nanoscope.png'
cell_picking_metadata_wildcard        = '{sample}/{modality}_{barcode}/cell_picking/metadata.csv'
matrix_bins_wildcard                  = '{sample}/{modality}_{barcode}/matrix/matrix_bin_{bins}/matrix.mtx.gz'
matrix_genebody_promoter_wildcard     = '{sample}/{modality}_{barcode}/matrix/matrix_genebody_promoter/matrix.mtx.gz'
fragments_noLA_wildcard               = '{sample}/{modality}_{barcode}/cellranger/outs/fragments_noLA_duplicates.tsv.gz'

# Bulk wildcards
macs_merged_per_modality_wildcard     = '{sample}/{modality}/peaks/macs2/{modality}_peaks_3column.bed'

rule all_preprocess:
    input:
    # Single-cell outputs
        cellranger           = [cellranger_bam_wildcard.format(sample=config['name'], modality=barcodes_dict[barcode],barcode=barcode) for barcode in barcodes_dict.keys()],
        cellranger_clean_out = [cellranger_cleanup_wildcard.format(sample=config['name'], modality=barcodes_dict[barcode],barcode=barcode) for barcode in barcodes_dict.keys()],
        cell_picking_out     = [cell_picking_metadata_wildcard.format(sample=config['name'], modality=barcodes_dict[barcode],barcode=barcode) for barcode in barcodes_dict.keys()],
        matrix_bins_out      = [matrix_bins_wildcard.format(sample=config['name'], modality=barcodes_dict[barcode],barcode=barcode, bins=bin_size) for barcode in barcodes_dict.keys() for bin_size in config['general']['bin_sizes']],
        genebody_matrix      = [matrix_genebody_promoter_wildcard.format(sample=config['name'], modality=barcodes_dict[barcode],barcode=barcode) for barcode in barcodes_dict.keys()],
        fragments_noLA       = [fragments_noLA_wildcard.format(sample=config['name'], modality=barcodes_dict[barcode],barcode=barcode) for barcode in barcodes_dict.keys()],
    # Bulk outputs
        bigwig       = expand('{sample}/{modality}/mapping_out/{modality}_merged.bw',sample = config['name'], modality = list(set(barcodes_dict.values()))),
        macs_3column = expand('{sample}/{modality}/peaks/macs2/{modality}_peaks.broadPeak',sample = config['name'], modality = list(set(barcodes_dict.values())))


rule create_barcodes_table:
    input:
        script = workflow.basedir + '/scripts/config2barcodes_table.py',
    output:
        barcodes_table = '{sample}/general/{sample}_barcodes_table.txt'
    params:
        sample_name = config['name'],
        config_path = workflow.configfiles[0], # assuming 1 configfile
    conda: '../envs/nanoscope_general.yaml'
    threads: 1
    shell:
        "python {input.script} --configfile {params.config_path} --sample {params.sample_name} --outfile {output.barcodes_table} --reversecomplement 2>&1"


rule split_R2:
    input:
        R2 = lambda wildcards: find_fastq_in_folder(path=config['fastq_path'], filters = [wildcards.lane,'R2']),
    output:
        modality = '{sample}/fastq_split_R2/{prefix}_{number}_{lane}_R2_001_modality.fastq.gz',
        singlecell = '{sample}/fastq_split_R2/{prefix}_{number}_{lane}_R2_001_singlecell.fastq.gz',
    params:
        outdir = '{sample}/fastq_split_R2/',
        script = workflow.basedir + '/scripts/process_R2_barcodes.py',
    conda: "../envs/nanoscope_general.yaml"
    threads: 1
    shell:
        "mkdir -p {params.outdir} && "
        "python3 {params.script} --input {input.R2} --output_folder {params.outdir}"

rule demultiplex:
    input:
        R1            = lambda wildcards: find_fastq_in_folder(path=config['fastq_path'], filters = [wildcards.lane,'R1']),
        R2_modality   = '{sample}/fastq_split_R2/{prefix}_{number}_{lane}_R2_001_modality.fastq.gz',
        R2_singlecell = '{sample}/fastq_split_R2/{prefix}_{number}_{lane}_R2_001_singlecell.fastq.gz',
        R3            = lambda wildcards: find_fastq_in_folder(path=config['fastq_path'], filters = [wildcards.lane, 'R3']),
        barcodes_table = '{sample}/general/{sample}_barcodes_table.txt',
    output:
        output = expand('{{sample}}/fastq_debarcoded/barcode_{barcode}/{{prefix}}_{{number}}_{{lane}}_{read}_{{suffix}}', barcode = config['barcodes'].keys(), 
                                                                                                                          read = ['R1','R2','R3']),
    params:
        outdir   = '{sample}/fastq_demultiplexed/',
        log    = "{sample}/fastq_debarcoded/{sample}_{lane}_log.txt",
        all_barcodes = ' '.join(config['barcodes'].keys()),
    conda: "../envs/nanoscope_general.yaml"
    threads: 8
    shell:
        """
        # Create folders for each barcode
        for b in {params.all_barcodes}; do
            mkdir -p {wildcards.sample}/fastq_debarcoded/barcode_$b;
        done
        # And for unmatched reads
        mkdir -p {wildcards.sample}/fastq_debarcoded/barcode_unmatched;
        # Demultiplex using fastq-multx
        fastq-multx {input.barcodes_table} {input.R2_modality} {input.R1} {input.R2_singlecell} {input.R3} \
        -o \
        {wildcards.sample}/fastq_debarcoded/barcode_%/{wildcards.prefix}_{wildcards.number}_{wildcards.lane}_modality_{wildcards.suffix} \
        {wildcards.sample}/fastq_debarcoded/barcode_%/{wildcards.prefix}_{wildcards.number}_{wildcards.lane}_R1_{wildcards.suffix} \
        {wildcards.sample}/fastq_debarcoded/barcode_%/{wildcards.prefix}_{wildcards.number}_{wildcards.lane}_R2_{wildcards.suffix} \
        {wildcards.sample}/fastq_debarcoded/barcode_%/{wildcards.prefix}_{wildcards.number}_{wildcards.lane}_R3_{wildcards.suffix} \
        > {params.log};
        """

rule run_cellranger:
    input:
        reads_list = lambda wildcards: get_fastq_for_cellranger(fastq_folder = config['fastq_path'], sample=config['name'], modality=wildcards.modality, barcode=wildcards.barcode),
        peaks      = macs_merged_per_modality_wildcard,
    output:
        bam   = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        frag  = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        meta  = '{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv',
        # peaks = '{sample}/{modality}_{barcode}/cellranger/outs/peaks.bed',
    params:
        cellranger_software = config['general']['cellranger_software'],
        cellranger_ref      = config['general']['cellranger_ref'],
        fastq_folder        = lambda wildcards: os.getcwd() + '/{sample}/fastq_debarcoded/barcode_{barcode}/'.format(sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode),
        sample_names        = lambda wildcards: ','.join(get_sample_names_for_cellranger(fastq_folder = config['fastq_path'])),
    threads: 8
    resources:
        mem_gb = 32
    shell:
        'rm -rf {wildcards.sample}/{wildcards.modality}_{wildcards.barcode}/cellranger/; '
        'cd {wildcards.sample}/{wildcards.modality}_{wildcards.barcode}/; '
        '{params.cellranger_software} count --id cellranger --reference {params.cellranger_ref} --fastqs {params.fastq_folder} --sample {params.sample_names} --peaks ../../{input.peaks} --localcores {threads} --localmem {resources.mem_gb}  2>&1 ; '

rule clean_cellranger_output:
    input:
        bam   = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        frag  = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        meta  = '{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv',
    output:
        '{sample}/{modality}_{barcode}/cellranger/outs/_cellranger_cleanup'
    params:
        cellranger_folder = str(Path('{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam').parents[1].resolve())
    shell:
        'touch {params.cellranger_folder}/tmp.txt;'                     # Create temp empty file to avoid error if the directory is empty
        'ls -d  {params.cellranger_folder}/* | grep -v outs | xargs rm -r; '
        'touch {output}'

rule barcode_metrics_peaks:
    input:
        bam   = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        peaks = '{sample}/{modality}/peaks/macs2/{modality}_peaks' + '_3column.bed'
    output:
        overlap_file = '{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt'
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
        peaks     = macs_merged_per_modality_wildcard,
        metadata  = cellranger_metadata_wildcard,
        fragments = cellranger_fragments_wildcard
    output:
        cell_picking_cells_10x_wildcard,
        cell_picking_cells_nanoscope_wildcard,
        cell_picking_metadata_wildcard,
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
    threads: 8
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
    threads: 8
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
    threads: 8
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

rule get_cells:
    input:
        metadata = '{sample}/{modality}_{barcode}/cell_picking/metadata.csv',
    output:
        cells = temp('{sample}/{modality}_{barcode}/matrix/all_cells.txt')
    params:
        script = workflow.basedir + '/scripts/get_passed_cells_barcodes.awk'
    shell:
        "cat {input} | awk -f {params.script} | sed 's/\"//g'| sort | uniq > {output.cells}"
        
rule create_matrix_bins:
    input:
        bam   = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        frag  = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        cells = '{sample}/{modality}_{barcode}/matrix/all_cells.txt',
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

rule download_annotation_and_get_genebody_and_promoter_gtf:
    output:
        gencode_gtf  = 'annotation/gencode.full.gtf.gz',
        genebody_gtf = 'annotation/genebody_and_promoter.gtf',
        gene_names   = 'annotation/gene_names.txt'
    params:
        url = config['general']['gtf'],
        script = workflow.basedir + '/scripts/get_filtered_genebody_promoter_from_gencode_gtf.awk',
        script2 = workflow.basedir + '/scripts/get_gene_name_from_gencode_gtf.awk'
    shell:
        "wget {params.url} -O {output.gencode_gtf}; "
        " gunzip -c {output.gencode_gtf} | awk -f {params.script} | cut -f1,4,5  > {output.genebody_gtf}; "
        " gunzip -c {output.gencode_gtf} | awk -f {params.script} | awk -f {params.script2} | sed 's/;//g' | sed 's/\"//g' > {output.gene_names} "

rule create_genebody_and_promoter_matrix:
    input:
        gtf   = 'annotation/genebody_and_promoter.gtf',
        frag  = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        cells = '{sample}/{modality}_{barcode}/matrix/all_cells.txt',
    output:
        features   = '{sample}/{modality}_{barcode}/matrix/matrix_genebody_promoter/features.tsv.gz',
        matrix     = '{sample}/{modality}_{barcode}/matrix/matrix_genebody_promoter/matrix.mtx.gz',
        barcodes   = '{sample}/{modality}_{barcode}/matrix/matrix_genebody_promoter/barcodes.tsv',
        folder     = directory('{sample}/{modality}_{barcode}/matrix/matrix_genebody_promoter/'),
    conda: '../envs/nanoscope_general.yaml'
    shell:
        'fragtk matrix -f {input.frag} -b {input.gtf} -c {input.cells} -o {output.folder}; '
        ''    
        




# rule bam_to_bw: # For QC reasons
#     input:
#         cellranger_bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'
#     output:
#         bigwig='{sample}/{modality}_{barcode}/bigwig/all_reads.bw'
#     threads: 8
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

