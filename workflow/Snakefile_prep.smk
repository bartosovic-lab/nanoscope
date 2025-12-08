import itertools
import os
import glob
import sys

# configfile: workflow.basedir + '/../config/config.yaml'

# samples_list  = list(config['samples'].keys())
# print(samples_list)

# barcodes_dict = {sample: config['samples'][sample]['barcodes'] for sample in samples_list}
# print(barcodes_dict)

# antibodies_list = list(set(itertools.chain(*[barcodes_dict[sample].keys() for sample in barcodes_dict.keys()])))
# print(antibodies_list)

# Find combinations of modalities
# These are only realistic combinations of modalities
# modalities_combinations =  [itertools.combinations(list(barcodes_dict[sample].keys()),i) for sample in samples_list for i in range(2,1+len(barcodes_dict[sample].keys()))]
# modalities_combinations = [list(x) for x in list(set(itertools.chain(*modalities_combinations)))]
# print(modalities_combinations)

# Define features
features = ['peaks']
bins = [5000,10000]

print(features)



def revcomp(seq: str) -> str:
    """Return reverse complement of an A/C/G/T sequence."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def get_fastq_for_cellranger(fastq_folder,sample,modality,barcode):
    import glob
    result = []
    fastq_folder = fastq_folder + "/**/*.fastq.gz"
    
    # Find all fastq files in a folder
    sys.stderr.write('Looking for fastq files in folder: {}\n'.format(fastq_folder))
    all_fastq_files  = glob.glob(fastq_folder,recursive=True)
    
    # Check if there are any fastq files and if they exist 
    check_fastq(all_fastq_files)
    
    # Parse
    all_fastq_parsed = [parse_fastq(x) for x in all_fastq_files]
    sys.stderr.write('Found {} fastq files: {}\n'.format(len(all_fastq_files),'\n'.join(all_fastq_files)))
    for x in all_fastq_parsed:
        if x['read'] == 'I1':
            continue
        result.append('{sample}/debarcoded_fastq/barcode_{barcode}/{sample}_{number}_{lane}_{read}_{suffix}'.format(\
            sample=sample, modality=modality , barcode=barcode, seq_id=x['id'], number=x['number'], lane=x['lane'], suffix=x['suffix'], read = x['read']))
    return(result)

def check_fastq(all_fastq_files):
    if len(all_fastq_files) == 0:
        sys.stderr.write('*** Error: Found 0 files in folder {}\n'.format(fastq_folder))
        sys.stderr.write('*** Aborting now! \n')
        raise Exception('No files found in fastq folder\n')
    for x in all_fastq_files:
        if not os.path.isfile(x):
            sys.stderr.write("*** Error: File {} does not exist\n".format(x))
            sys.stderr.write('*** Aborting now! \n')
            raise Exception("File does not exist\n")
    return

def parse_fastq(path):
    import os
    import re
    result = {}
    fastq = os.path.basename(path)
    result['number'] = re.findall('_S[0-9]+_', fastq)[0].strip("_")
    result['lane']   = re.findall('_L[0-9]+_', fastq)[0].strip("_")
    result['read']   = re.findall('_[RI][0-9]+_', fastq)[0].strip("_")
    result['id']     = re.split('_S[0-9]+_',fastq)[0].strip("_")
    result['suffix']  = re.split('_[RI][0-9]+_',fastq)[1].strip("_")
    return(result)

# Check if program is installed

def check_installed(path):
    import os
    if os.system("command -v " + path) == 0:
        return True
    else:
        return False


def generate_matrix_out(sample, modality, barcode, bins = bins):
    if not check_installed('fragtk'):
        sys.stderr.write('*** Warning: fragtk not installed, the pipeline will not generate matrix files\n')
        return []
    
    matrix_peaks = ['{sample}/{modality}_{barcode}/matrix/matrix_peaks/'.format(sample=sample, modality=modality,barcode=barcode)]
    matrix_bins  = ['{sample}/{modality}_{barcode}/matrix/matrix_bin_{bins}/'.format(sample=sample, modality=modality,barcode=barcode, bins = b) for b in bins]
    matrix_genes = ['{sample}/{modality}_{barcode}/matrix/matrix_genes/'.format(sample=sample, modality=modality,barcode=barcode)]
    print(matrix_genes)
    return matrix_peaks + matrix_bins + matrix_genes
