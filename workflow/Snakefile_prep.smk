import itertools
import os
import glob
import re
import collections
from pathlib import Path

# Some variable
bwa_index                          = str(Path(config['general']['cellranger_ref'] + '/fasta/genome.fa'))

########## All of the wildcards used in the pipeline ##########
# Bulk wildcards
debarcoded_fastq_wildcard               = '{sample}/fastq_debarcoded/barcode_{barcode}/{prefix}_{number}_{lane}_{read}_{suffix}'
trimmed_fastq_wildcard                  = '{sample}/{modality}_{barcode}/fastq_trimmed/{prefix}_{number}_{lane}_{read}_{suffix}'
bowtie2_map_wildcard                    = '{sample}/{modality}_{barcode}/mapping_out/{sample}_{modality}_{lane}_mapped.bam'
bwa_map_wildcard                        = '{sample}/{modality}_{barcode}/mapping_out/{sample}_{modality}_{lane}_mapped.bam'
bam_sorted_wildcard                     = '{sample}/{modality}_{barcode}/mapping_out/{sample}_{modality}_{lane}_sorted.bam'
macs_wildcard                           = '{sample}/{modality}_{barcode}/peaks/macs2/{sample}_{modality}_peaks.broadPeak'
bam_merged_wildcard                     = '{sample}/{modality}/mapping_out/{modality}_merged.bam'
bigwig_wildcard                         = '{sample}/{modality}/mapping_out/{modality}_merged.bw'
macs_merged_per_modality_wildcard       = '{sample}/{modality}/peaks/macs2/{modality}_peaks.broadPeak'
fasta_index_wildcard                    = 'fasta_index.fai'
# bowtie2_index_wildcard             = '{sample}/reference/bowtie2/genome'

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

def find_all_fastq_files(path):
    all_fastq = itertools.chain(*[glob.glob(path + x) for x in ['/**/*R*.fastq.gz', '/*R*.fastq.gz']])
    all_fastq = [fastq_file(x) for x in all_fastq]
    return (all_fastq)

def invert_dict(d):
    d_new = {}
    for k,v in d.items():
        if v not in d_new:
            d_new[v] = [k]
        else:
            d_new[v].append(k)
    return(d_new)

def mask_wildcards(wildcard_string, list_of_wildcards, invert=False):
    all_wildcards = re.findall('{(.*?)}',wildcard_string)
    for wildcard in all_wildcards:
        # Forward masking - only the ones in the list
        if wildcard in list_of_wildcards and not invert:
            wildcard_string = wildcard_string.replace('{' + wildcard + '}', '{{' + wildcard + '}}')
        # Invert masking - all but the ones in the list
        elif wildcard not in list_of_wildcards and invert:
            wildcard_string = wildcard_string.replace('{' + wildcard + '}', '{{' + wildcard + '}}')            
    return wildcard_string

class snakemake_run:
    def __init__(self,config):
        self.tempdir = config['general']['tempdir']
        self.cellranger_software = config['general']['cellranger_software']
        self.cellranger_ref = config['general']['cellranger_ref']
        self.macs_genome = config['general']['macs_genome']
        self.samples_list = list(config['samples'].keys())
        self.features = 'peaks'
        self.samples = {}
        for s in self.samples_list:
            self.samples[s] = sample(config,s)

class sample:
    def __init__(self,config,sample_name):
        self.sample_name     = sample_name
        self.fastq_path      = config['samples'][self.sample_name]['fastq_path']
        # barcodes_dict = { barcode1: modality1, barcode2: modality1, ...}
        self.barcodes_dict   = config['samples'][self.sample_name]['barcodes']
        # reverese_barcodes_dict = {modality1: [barcode1, barcode2, ...], modality2: [barcode1, barcode2, ...]}
        self.reverse_barcodes_dict = invert_dict(self.barcodes_dict)
        # modality_names = [modality1, modality2, ...]
        self.modality_names  = list(set([x for x in config['samples'][self.sample_name]['barcodes'].values()]))       # List of modalities
        # barcodes_list = [barcode1, barcode2, ...]
        self.barcodes_list   = list(set([x for x in config['samples'][self.sample_name]['barcodes'].keys()]))     # List of barcodes
        self.bins = ['5000','10000','25000','50000','100000']

        # Some checks here
        self.check_valid_barcodes()

        self.all_fastq_files = find_all_fastq_files(self.fastq_path)
        self.all_lanes       = sorted(list(set([x.lane for x in self.all_fastq_files])))
        self.all_reads       = sorted(list(set([x.read for x in self.all_fastq_files])))

        # Easy access to lane and read information through dictionary
        self.fastq_by_lane = collections.defaultdict(dict)
        for f in self.all_fastq_files:
            self.fastq_by_lane[f.lane][f.read] = f

        ###########
        # OUTPUTS #
        ###########
        # All outputs are based on wildcards declared above and are generated here

        self.generate_debarcoded_output(files_list='debarcoded_fastq_all', files_dict='debarcoded_fastq_dict', files_by_modality = 'debarcoded_fastq_by_modality', wildcard=debarcoded_fastq_wildcard)
        self.generate_debarcoded_output(files_list='trimmed_fastq_all', files_dict='trimmed_fastq_dict', files_by_modality = 'trimmed_fastq_by_modality', wildcard=trimmed_fastq_wildcard,filter_read = 'R2')
        
        # Creates following: 
        # debarcoded_fastq_all = [fastq_debarcoded_1.fastq, fastq_debarcoded_2.fastq, ...]
        # debracoded_fastq_dict = {lane: {barcode: {read: fastq_file1, fastq_file2, ...}, barcode2: {read: fastq_file1, fastq_file2, ...}}, lane2: ...}
        # debarcoded_fastq_by_modality = {modality1: [fastq_file1, fastq_file2, ...], modality2: [fastq_file1, fastq_file2, ...]}


        # Debarcoded fastq files by modality

        # Bulk outputs
        self.bowtie2_bam_all                    = [bowtie2_map_wildcard.format(sample = self.sample_name,modality=self.barcodes_dict[b],barcode=b,lane=l) for b in self.barcodes_list for l in self.all_lanes]
        self.bam_sorted_all                     = [bam_sorted_wildcard.format(sample = self.sample_name,modality=self.barcodes_dict[b],barcode=b,lane=l) for b in self.barcodes_list for l in self.all_lanes]
        self.bam_merged_all                     = [bam_merged_wildcard.format(sample=self.sample_name,modality=self.barcodes_dict[b],barcode=b) for b in self.barcodes_list]
        self.bigwig_all                         = [bigwig_wildcard.format(sample=self.sample_name,modality=self.barcodes_dict[b],barcode=b) for b in self.barcodes_list]
        self.macs_all                           = [macs_wildcard.format(sample=self.sample_name,modality=self.barcodes_dict[b],barcode=b) for b in self.barcodes_list]
        self.macs_merged_accross_modality_all   = [macs_merged_per_modality_wildcard.format(sample=self.sample_name, modality = m) for m in self.modality_names]


        # Single-cell outputs
        self.cellranger_all         = [cellranger_fragments_wildcard.format(sample=self.sample_name,modality=self.barcodes_dict[b],barcode=b)  for b in self.barcodes_list] + \
                                      [cellranger_bam_wildcard.format(sample=self.sample_name,modality=self.barcodes_dict[b],barcode=b)  for b in self.barcodes_list]
        self.cellranger_cleanup_all = [cellranger_cleanup_wildcard.format(sample=self.sample_name,modality=self.barcodes_dict[b],barcode=b)  for b in self.barcodes_list]
        self.cell_picking_all       = [cell_picking_metadata_wildcard.format(sample=self.sample_name,modality=self.barcodes_dict[b],barcode=b)  for b in self.barcodes_list]
        self.matrix_bins_all        = [matrix_bins_wildcard.format(sample=self.sample_name,modality=self.barcodes_dict[b],barcode=b,bins=bin)  for b in self.barcodes_list for bin in self.bins]
        self.matrix_genebody_all    = [matrix_genebody_promoter_wildcard.format(sample=self.sample_name,modality=self.barcodes_dict[b],barcode=b)  for b in self.barcodes_list]

    # Takes a list of files and returns dictionary of files dict[lane][barcode][read]
    def generate_debarcoded_output(self, files_list, files_dict, files_by_modality, wildcard,filter_read = False):
        setattr(self, files_list,[])    # Empty list
        setattr(self,files_dict, {l: collections.defaultdict(dict) for l in self.all_lanes})    # Empty dictionary

        for f in self.all_fastq_files:                  # For each file in raw input
            for barcode in self.barcodes_list:        # Multiply for each file by barcode
                if f.read == filter_read:
                    continue
                fastq_path     = wildcard.format(sample = self.sample_name, modality = self.barcodes_dict[barcode], barcode = barcode, prefix = f.prefix, number = f.number, lane = f.lane, read = f.read, suffix = f.suffix)

                # Dictionary to access the files by lane, modality and read
                d = getattr(self,files_dict)
                d[f.lane][barcode][f.read] = fastq_file(fastq_path, modality = self.barcodes_dict[barcode], barcode = barcode)
                setattr(self,files_dict, d)

                # List of all files generated
                d=getattr(self,files_list)
                d.append(getattr(self,files_dict)[f.lane][barcode][f.read])
                setattr(self,files_list,d)
                setattr(self, files_by_modality, {m: [x for x in getattr(self,files_list) if x.modality == m] for m in self.modality_names})

        return(self)

    def check_valid_barcodes(self, alphabet = 'ATCG'):
        for b in self.barcodes_list:
            if not set(b).issubset(set(alphabet)):
                raise ValueError('Barcode {b} contains non-ATCG characters\n{barcodes}\n'.format(b=b,barcodes = '\n'.join(self.barcodes_list)))

class fastq_file:
    def __init__(self,path,modality = False,barcode=False):
        self.path     = path
        self.abspath  = os.path.abspath(path)
        self.basename = os.path.basename(self.abspath)
        self.dirname  = os.path.dirname(self.abspath)
        self.modality = modality
        self.barcode  = barcode
        self.parse_fastq(self.basename)

    def parse_fastq(self, path):

        # Parse fastq file name into following [e.g. P29308_1002_S2_L002_R3_001.fastq.gz]
        # 1. fastq_prefix   [P29308_1002]
        # 2. fastq_number   [S2]
        # 3. lane number    [L002]
        # 4. read number    [R3]
        # 5. fastq suffix   [_001.fastq.gz]

        self.prefix = re.split('_S[0-9]+_',path)[0].strip("_")
        self.number = re.findall('_S[0-9]+_',path)[0].strip("_")
        self.lane   = re.findall('_L[0-9]+_',path)[0].strip("_")
        self.read   = re.findall('_[RI][0-9]+_',path)[0].strip("_")
        self.suffix = re.split('_[RI][0-9]+_',path)[1].strip("_").replace('.gz','')

        return(self)
    
    def __repr__(self):
        return(self.path)

run = snakemake_run(config)

debarcoded_fastq_output = [debarcoded_fastq_wildcard.format(sample = '{sample}', barcode = b, prefix = '{prefix}', number = '{number}', lane = '{lane}', read = r, suffix = '{suffix}') for s in run.samples_list for b in run.samples[s].barcodes_dict.keys() for r in ['R1','R2','R3']]
trimmed_fastq_output    = {r: trimmed_fastq_wildcard.replace('{read}',r) for r in ['R1','R2','R3']}

