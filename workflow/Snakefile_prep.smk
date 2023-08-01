import itertools
import os
import glob
import re
import collections

debarcoded_fastq_wildcard = '{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/{prefix}_{number}_{lane}_{read}_{suffix}'
trimmed_fastq_wildcard    = {'R1': '{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/trimmed/{prefix}_{lane}_val_R1.fq.gz',
                             'R2': '{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/trimmed/{prefix}_{lane}_val_R2.fq.gz'}

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
        print(self.all_reads)


        # Easy access to lane and read information through dictionary
        self.fastq_by_lane = collections.defaultdict(dict)
        for f in self.all_fastq_files:
            self.fastq_by_lane[f.lane][f.read] = f

        #########################################
        # Files that are output of debarcode.py #
        #########################################
        self.generate_debarcode_output(debarcoded_fastq_wildcard)

    def generate_debarcode_output(self,debarcoded_fastq_wildcard):
        self.debarcoded_fastq_by_lane = {}
        self.debarcoded_fastq_all = []
        for f in self.all_fastq_files:
            self.debarcoded_fastq_by_lane[f.lane] = collections.defaultdict(dict)
            for modality in self.modality_names:
                fastq_path     = debarcoded_fastq_wildcard.format(sample = self.sample_name, modality = modality, barcode = self.barcodes_dict[modality], prefix = f.prefix, number = f.number, lane = f.lane, read = f.read, suffix = f.suffix)
                self.debarcoded_fastq_by_lane_wildcard                  = debarcoded_fastq_wildcard
                self.debarcoded_fastq_by_lane[f.lane][modality][f.read] = fastq_file(fastq_path)
                self.debarcoded_fastq_all.append(fastq_file(fastq_path))
        return(self)
    #
    # def generate_trim_output(self):
    #     self.all_trimmed_fastq = []
    #     self.trimmed_fastq_by_lane ={}
    #     for f in self.debarcoded_fastq_all:
    #
    #
    #
    #




class fastq_file:
    def __init__(self,path):
        self.path     = path
        self.abspath  = os.path.abspath(path)
        self.basename = os.path.basename(self.abspath)
        self.dirname  = os.path.dirname(self.abspath)
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



# samples_list  = list(config['samples'].keys())
# print(samples_list)
#
# barcodes_dict = {sample: config['samples'][sample]['barcodes'] for sample in samples_list}
# print(barcodes_dict)
#
# # Define features
# features = ['peaks']
# print(features)

# print("Class attr from here")
# print(run.samples_list)
# print(run.features)
# print(run.samples['BARNYARD_1'].lanes)

# Get lanes in this experiment (REGEX MATCH _L[0-9]+_ and remove underscores)
# lanes = [glob.glob(config['samples'][sample]['fastq_path'] + x) for x in ['/**/*R1*.fastq.gz','/*R1*.fastq.gz'] for sample in config['samples'].keys()]
# lanes = [os.path.basename(x) for x in itertools.chain(*lanes)]
# lanes = sorted(list(set([x[re.search('_L[0-9]+_',x).start()+1:re.search('_L[0-9]+_',x).end()-1] for x in lanes])))
# print(lanes)


#
# def get_fastq_for_cellranger(fastq_folder,sample,modality,barcode):
#     import glob
#     result = []
#     all_fastq_files  = glob.glob(fastq_folder + "/**/*.fastq.gz",recursive=True)
#     all_fastq_parsed = [parse_fastq(x) for x in all_fastq_files]
#     for x in all_fastq_parsed:
#         if x['read'] == 'I1':
#             continue
#         result.append('{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/{sample}_{number}_{lane}_{read}_{suffix}'.format(\
#             sample=sample, modality=modality , barcode=barcode, seq_id=x['id'], number=x['number'], lane=x['lane'], suffix=x['suffix'], read = x['read']))
#     return(result)
