import itertools
import os
import glob
import re
import collections
from pathlib import Path

debarcoded_fastq_wildcard = '{sample}/{modality}_{barcode}/fastq_debarcoded/barcode_{barcode}/{prefix}_{number}_{lane}_{read}_{suffix}'
trimmed_fastq_wildcard    = '{sample}/{modality}_{barcode}/fastq_trimmed/{prefix}_{number}_{lane}_{read}_{suffix}'
bowtie2_map_wildcard      = '{sample}/{modality}_{barcode}/bowtie2_out/{sample}_{modality}_{lane}_mapped.bam'
bam_sorted_wildcard       = '{sample}/{modality}_{barcode}/bowtie2_out/{sample}_{modality}_{lane}_sorted.bam'
bam_merged_wildcard       = '{sample}/{modality}_{barcode}/bowtie2_out/{sample}_{modality}_merged.bam'
bigwig_wildcard           = '{sample}/{modality}_{barcode}/bowtie2_out/{sample}_{modality}_merged.bw'

debarcoded_fastq_output = {r: '{sample}/{modality}_{barcode}/fastq_debarcoded/barcode_{barcode}/{prefix}_{number}_{lane}_{read}_{suffix}'.replace('{read}',r) for r in ['R1','R2','R3']}
trimmed_fastq_output    = {r: '{sample}/{modality}_{barcode}/fastq_trimmed/{prefix}_{number}_{lane}_{read}_{suffix}'.replace('{read}',r) for r in ['R1','R2','R3']}

debarcode_params_outdir = str(Path(debarcoded_fastq_wildcard).parents[1])
trim_params_outdir      = str(Path(trimmed_fastq_wildcard).parents[0])

def find_all_fastq_files(path):
    all_fastq = itertools.chain(*[glob.glob(path + x) for x in ['/**/*R*.fastq.gz', '/*R*.fastq.gz']])
    all_fastq = [fastq_file(x) for x in all_fastq]
    return (all_fastq)


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
        self.barcodes_dict   = config['samples'][self.sample_name]['barcodes']
        self.modality_names  = [x for x in config['samples'][self.sample_name]['barcodes'].keys()]
        self.barcodes_list   = [x for x in config['samples'][self.sample_name]['barcodes'].values()]

        self.all_fastq_files = find_all_fastq_files(self.fastq_path)
        self.all_lanes       = sorted(list(set([x.lane for x in self.all_fastq_files])))
        self.all_reads       = sorted(list(set([x.read for x in self.all_fastq_files])))
        # print(self.all_reads)


        # Easy access to lane and read information through dictionary
        self.fastq_by_lane = collections.defaultdict(dict)
        for f in self.all_fastq_files:
            self.fastq_by_lane[f.lane][f.read] = f

        ###########
        # OUTPUTS #
        ###########

        self.generate_debarcoded_output(files_list='debarcoded_fastq_all', files_dict='debarcoded_fastq_dict', wildcard=debarcoded_fastq_wildcard)
        self.generate_debarcoded_output(files_list='trimmed_fastq_all', files_dict='trimmed_fastq_dict',wildcard=trimmed_fastq_wildcard,filter_read = 'R2')

        # Bowtie2 mapping output
        self.bowtie2_bam_all  = [bowtie2_map_wildcard.format(sample = self.sample_name,modality=m,barcode=self.barcodes_dict[m],lane=l) for m in self.modality_names for l in self.all_lanes]
        self.bam_sorted_all   = [bam_sorted_wildcard.format(sample = self.sample_name,modality=m,barcode=self.barcodes_dict[m],lane=l) for m in self.modality_names for l in self.all_lanes]
        self.bam_merged_all   = [bam_merged_wildcard.format(sample=self.sample_name,modality=m,barcode=self.barcodes_dict[m]) for m in self.modality_names]
        self.bigwig_all       = [bigwig_wildcard.format(sample=self.sample_name,modality=m,barcode=self.barcodes_dict[m]) for m in self.modality_names]


    def generate_debarcoded_output(self, files_list, files_dict, wildcard,filter_read = False):
        setattr(self, files_list,[])    # Empty list
        setattr(self,files_dict, {l: collections.defaultdict(dict) for l in self.all_lanes})    # Empty dictionary

        for f in self.all_fastq_files:                  # For each file in raw input
            for modality in self.modality_names:        # Multiply for each file by modality
                if f.read == filter_read:
                    continue
                fastq_path     = wildcard.format(sample = self.sample_name, modality = modality, barcode = self.barcodes_dict[modality], prefix = f.prefix, number = f.number, lane = f.lane, read = f.read, suffix = f.suffix)

                # Dictionary to access the files by lane, modality and read
                d = getattr(self,files_dict)
                d[f.lane][modality][f.read] = fastq_file(fastq_path, modality = modality)
                setattr(self,files_dict, d)

                # List of all files generated
                d=getattr(self,files_list)
                d.append(getattr(self,files_dict)[f.lane][modality][f.read])
                setattr(self,files_list,d)
        return(self)

class fastq_file:
    def __init__(self,path,modality = False):
        self.path     = path
        self.abspath  = os.path.abspath(path)
        self.basename = os.path.basename(self.abspath)
        self.dirname  = os.path.dirname(self.abspath)
        self.modality = modality
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
        self.suffix = re.split('_[RI][0-9]+_',path)[1].strip("_")

        return(self)

run = snakemake_run(config)