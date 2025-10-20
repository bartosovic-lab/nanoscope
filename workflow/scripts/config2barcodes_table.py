import yaml
import sys
import argparse

parser = argparse.ArgumentParser(description='Generate barcodes table from config file')
parser.add_argument('--configfile', type=str, help='Path to the config YAML file')
parser.add_argument('--sample', type=str, help='Sample name')
parser.add_argument('--reversecomplement', action='store_true', help='Whether to reverse complement the barcodes')
parser.add_argument('--outfile', type=str, help='Output file path')
args = parser.parse_args()

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

with open(args.configfile, 'r') as file:
    config = yaml.safe_load(file)
try: 
    barcodes_table = config['barcodes']
except KeyError:
    print(f"Error: Sample '{args.sample}' not found in the configuration file.")
    sys.exit(1)

with open(args.outfile, 'w') as out:
    print(barcodes_table)
    for key in barcodes_table:
        sample = key
        if args.reversecomplement:
            key = reverse_complement(key)
        out.write("{s}\t{b}\n".format(b=key, s=sample))
        #out.write("{s}\t{b}\t{s}_{b}_{lane}_R1.fastq.gz\t{s}_{b}_{lane}_R2.fastq.gz\n".format(b=key, s=barcodes_table[key], lane = lane))






