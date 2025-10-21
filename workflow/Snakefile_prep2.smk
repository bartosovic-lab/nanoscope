import itertools
import os
import glob
import sys

# configfile: workflow.basedir + '/../config/config.yaml'


barcodes_dict = config['barcodes']
print(barcodes_dict)

antibodies_list = list(set(itertools.chain(*[barcodes_dict.values()])))
print(antibodies_list)

# Find combinations of modalities
# These are only realistic combinations of modalities
#modalities_combinations =  [itertools.combinations(list(barcodes_dict[sample].keys()),i) for sample in samples_list for i in range(2,1+len(barcodes_dict[sample].keys()))]
#modalities_combinations = [list(x) for x in list(set(itertools.chain(*modalities_combinations)))]
#print(modalities_combinations)

# Define features
features = ['peaks']
bins = [5000,10000]

print(features)

def find_fastq_in_folder(path,filters = ['']):
    import glob
    all_fastq = glob.glob(path + "/**/*.fastq.gz", recursive=True)
    matched = [f for f in all_fastq if all(x in f for x in filters)]
    if len(matched) == 1:
        return matched[0]
    elif len(matched) > 1:
        raise ValueError(f"Multiple files matched for filters {filters}, found: {matched}, path: {path}")
    else:
        raise ValueError(f"No files matched for filters {filters}, found: {matched}, path: {path}")

def merge_bam_inputs(fastq_folder,sample,modality,barcodes_dict):
    import glob
    result = []
    fastq_folder = fastq_folder + "/**/*.fastq.gz"
    all_fastq_files = glob.glob(fastq_folder,recursive=True)
    check_fastq(all_fastq_files)
    all_fastq_parsed = [parse_fastq(x) for x in all_fastq_files]
    for x in all_fastq_parsed:
        if x['read'] != 'R1':
            continue
        for barcode in barcodes_dict:
            if barcodes_dict[barcode] != modality:
                continue
            result.append('{sample}/{modality}_{barcode}/mapping_out/{prefix}_{number}_{lane}_{suffix}_sorted.bam'.format(\
                sample=sample, modality=modality , barcode=barcode, prefix=x['id'], number=x['number'], lane=x['lane'], suffix=x['suffix']))
    return(result)


def parse_fastq(path):
    import os
    import re
    result = {}
    fastq = os.path.basename(path)
    result['id']     = re.split('_S[0-9]+_',fastq)[0].strip("_")
    result['number'] = re.findall('_S[0-9]+_', fastq)[0].strip("_")
    result['lane']   = re.findall('_L[0-9]+_', fastq)[0].strip("_")
    result['read']   = re.findall('_[RI][0-9]+_', fastq)[0].strip("_")
    result['suffix']  = re.split('_[RI][0-9]+_',fastq)[1].strip("_")
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

def get_fastq_for_cellranger(fastq_folder,sample,modality,barcode):
    # This retrieves filenames of fastq files, and reformats them to the cellranger expected format after demultiplexing
    # Also this returns a list of files, so that all parallel demultiplexing outputs are included in single cellranger run
    import glob
    result = []
    fastq_folder = fastq_folder + "/**/*.fastq.gz"
    
    # Find all fastq files in a folder
    # sys.stderr.write('Looking for fastq files in folder: {}\n'.format(fastq_folder))
    all_fastq_files = glob.glob(fastq_folder,recursive=True)
    
    # Check if there are any fastq files and if they exist 
    check_fastq(all_fastq_files)
    
    # Parse
    all_fastq_parsed = [parse_fastq(x) for x in all_fastq_files]
    # sys.stderr.write('Found {} fastq files: {}\n'.format(len(all_fastq_files),'\n'.join(all_fastq_files)))
    for x in all_fastq_parsed:
        if x['read'] == 'I1':
            continue
        result.append('{sample}/fastq_debarcoded/barcode_{barcode}/{seq_id}_{number}_{lane}_{read}_{suffix}'.format(\
            sample=sample, modality=modality , barcode=barcode, seq_id=x['id'], number=x['number'], lane=x['lane'], suffix=x['suffix'], read = x['read']))
    return(result)

def get_sample_names_for_cellranger(fastq_folder):
    import glob
    fastq_folder = fastq_folder + "/**/*.fastq.gz"
    all_fastq_files = glob.glob(fastq_folder,recursive=True)
    check_fastq(all_fastq_files)
    all_fastq_parsed = [parse_fastq(x) for x in all_fastq_files]
    sample_names = list(set([x['id'] for x in all_fastq_parsed]))
    if sample_names == []:
        sys.stderr.write('*** Error: No sample names found in fastq folder {}\n'.format(fastq_folder))
        sys.stderr.write('*** Aborting now! \n')
        raise Exception('No sample names found in fastq folder\n')
    return(sample_names)

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