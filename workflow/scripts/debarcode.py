import os
import argparse
from glob import glob
import sys
from re import split
import regex
import gzip
from contextlib import ExitStack
from collections import defaultdict
import time

import yaml
from pysam import FastxFile
import Levenshtein

def log(message):
    sys.stderr.write(time.strftime("%Y-%m-%d %H:%M:%S") + " " + message + "\n")

class bcdCT:
    def __init__(self,args):
        self.detect_input(args.input)
        self.detect_reads()
        self.single_cell=args.single_cell
        self.out_prefix=args.out_prefix
        if self.single_cell:
            self.out_reads = ['R1','R2','R3']
        else:
            self.out_reads = ['R1','R3']

        if args.name:
            self.name = args.name
        else:
            self.autodetect_name()


        self.autodetect_barcodes(args)
        self.prep_out_filenames()

    def detect_input(self,input):
        Error_message="*** Error: Wrong input files specified. The input must be either folder with _R1_*.fastq.gz _R2_*.fastq.gz _R3_*.fastq.gz files or paths to the files themselves ***" +\
        "The files should be placed in the same folder" +\
        "e.g. /data/path_to_my_files/*L001*.fastq.gz or /data/path_to_my_files/"

        input = [os.path.abspath(x) for x in input]
        if len(input) == 1 and os.path.isdir(input[0]):     # Case input is single directory
            self.input_dir = input[0]
            self.input_files = []
            self.input_files.extend(glob(self.input_dir + "/*.fastq.gz"))
            self.input_files.extend(glob(self.input_dir + "/*.fq.gz"))

        elif len(input) > 1:                                  # Case input are multiple files
            self.input_files = input
            self.input_dir = list(set([os.path.dirname(x) for x in self.input_files]))
            if not len(self.input_dir) == 1:
                log(Error_message)
                sys.exit(1)
            if not sum([x.endswith('.fastq.gz') or x.endswith('.fq.gz') for x in self.input_files]) == len(self.input_files):
                sys.exit(1)
                log(Error_message)
        else:
            sys.exit(1)
            log(Error_message)

    def detect_reads(self):
        Error_message="*** Error: Please specify exactly one _R1_ _R2_ and _R3_ file or folder with exactly one of each files ***" + \
                      "e.g. /data/path_to_my_files/*L001*.fastq.gz or /data/path_to_my_files/"
        self.path_in = {}
        self.path_in['R1'] = [x for x in self.input_files if "_R1_" in x]
        self.path_in['R2'] = [x for x in self.input_files if "_R2_" in x]
        self.path_in['R3'] = [x for x in self.input_files if "_R3_" in x]

        if len(self.path_in['R1']) != 1 or len(self.path_in['R2']) != 1 or len(self.path_in['R3']) != 1:
            log(Error_message)
            sys.exit(1)

        self.path_in = {key:self.path_in[key][0] for key in self.path_in.keys()}



    def in_handles(self,stack):
        in_stack = {x: stack.enter_context(FastxFile(self.path_in[x],'r')) for x in ['R1','R2','R3']}
        return in_stack

    def prep_out_filenames(self):
        self.path_out = {barcode: {} for barcode in self.picked_barcodes}
        self.path_out  = {barcode: {read: "{0}/barcode_{1}/{2}".format(self.out_prefix,barcode,os.path.basename(self.path_in[read])) for read in self.out_reads} for barcode in self.picked_barcodes}
        # If args.name is specified, replace the sample_id prefix with the one specified in args.name
        # e.g. nanoCT_MB22_001_S1_L001_R1_001.fastq.gz is input --name is test
        # Change to test_S1_L001_R1_001.fastq.gz
        if args.name:
            for barcode in self.path_out:
                for read in self.path_out[barcode]:
                    sample_id = split('_S[0-9]+_', os.path.basename(self.path_out[barcode][read]))[0].strip("_")
                    self.path_out[barcode][read] = self.path_out[barcode][read].replace(sample_id,args.name)

    def autodetect_name(self):
        Error_message = "*** Error: Prefix for R1 R2 R3 files not the same. Please use the same prefix for all the files or specify experiment name ***"

        self.name = [split("_R[0-9]_", str(x)) for x in self.path_in.values()]
        self.name = [x[0] for x in self.name]

        if len(list(set(self.name))) > 1:
            log(Error_message)
            sys.exit(1)

        self.name = self.name[0].split("/")[-1]

    def create_out_handles(self,stack):
        for bcd in self.picked_barcodes:
            os.makedirs(self.out_prefix + "/barcode_" + bcd, exist_ok=True)
        self.out_stack = {barcode: {read: stack.enter_context(gzip.open(self.path_out[barcode][read],'wt'))for read in self.out_reads} for barcode in self.picked_barcodes}


    def __iter__(self):
        with FastxFile(self.path_in['R1']) as f1, FastxFile(self.path_in['R2']) as f2, FastxFile(self.path_in['R3']) as f3:
            for r1,r2,r3 in zip(f1,f2,f3):
                yield r1, r2, r3

    def autodetect_barcodes(self,args):
        barcodes = defaultdict(int)
        n=0
        for read1,read2,read3 in self:
            hit     = find_seq(args.pattern, read2.sequence, nmismatch=0)
            MeA_hit = find_seq(pattern = args.no_barcode_seq, DNA_string=read2.sequence, nmismatch=2)
            if MeA_hit:
                barcodes['MeA'] += 1
            elif not hit or hit == 'Multiple':
                continue
            else:
                hit = int(hit)
                read_barcode = get_read_barcode(read2, hit)
                try:
                    barcodes[read_barcode] += 1
                except KeyError:
                    barcodes[read_barcode] = 1
            n += 1
            if n == 50000:
                break

        top_barcodes = sorted(barcodes, key=barcodes.get, reverse=True)[:args.Nbarcodes]
        picked_barcodes = {key: barcodes[key] for key in top_barcodes}
        log("Detected following most abundant barcodes out of first {} barcodes:\n{}".format(n, picked_barcodes))
        if args.barcode != "None":
            self.picked_barcodes = args.barcode
            log("Barcode specified for demultiplexing [{barcode}] in top found barcodes: {bool} ".format(bool = [(x,x in picked_barcodes.keys()) for x in args.barcode], barcode = args.barcode))
        else:
            self.picked_barcodes = [i for i in picked_barcodes.keys()]

        if args.report_MeA:
            self.picked_barcodes.append('MeA')
            log("MeA sequence will be reported in the output files due to --report_MeA flag")
        if args.report_no_hit:
            self.picked_barcodes.append('no_spacer')
            log("No spacer sequence will be reported in the output files due to --report_no_hit flag")
        
        print('final barcodes used for demultiplexing:')
        print(self.picked_barcodes)
        

def get_read_barcode(string,index):
    read_barcode = revcompl(string.sequence[index - 8:index])  # Get the barcode sequence
    return read_barcode

def extract_cell_barcode(read,index):
    read.sequence = read.sequence[index:index + 16]  # Get the cell barcode
    read.quality  = read.quality[index:index + 16]   # Get corresponding Quality score
    return read

def revcompl(seq):
    revcomp_table = {
        "A": "T",
        "G": "C",
        "C": "G",
        "T": "A",
        "N": "N"
    }
    complement = "".join([revcomp_table[letter] for letter in seq.upper()])  # Complement
    return complement[::-1]  # Reverse

def rev(seq):
    return seq[::-1]

def find_seq(pattern, DNA_string, nmismatch=2):
    for n in range(0,nmismatch + 1):
        r = regex.compile('({0}){{e<={1}}}'.format(pattern, n))
        res = r.finditer(DNA_string)
        hit = [x.start() for x in res]
        if len(hit) == 0:
            continue
        if len(hit) > 1:
            return None
        if len(hit) == 1:
            return int(hit[0])
    return None

def main(args):
    exp = bcdCT(args)
    statistics = defaultdict(int)
    log("Creating file output handles ")
    with ExitStack() as stack:
        exp.create_out_handles(stack)
        n = 0
        log("Starting demultiplexing ")
        for read1,read2,read3 in exp:
            n+=1
            if n % 5000000 == 0:
                log("{} reads processed".format(n))
            assert (read1.name == read2.name == read3.name)                                                 # Make sure the fastq files are ok

            spacer_hit = find_seq(pattern=args.pattern,DNA_string=read2.sequence,nmismatch=2)
            MeA_hit    = find_seq(pattern = args.no_barcode_seq, DNA_string=read2.sequence, nmismatch=2)
            
            if not spacer_hit and MeA_hit:
                read_barcode = get_read_barcode(read2, MeA_hit)                                               # Returns only barcode e.g. ACTGACTG
                hit_barcode  = 'MeA'
                if exp.single_cell:
                    read2 = extract_cell_barcode(read2, MeA_hit-16)     # The cell barcode is 16bp long and is positioned before the MeA spacer

            elif spacer_hit:
                read_barcode = get_read_barcode(read2, spacer_hit)                                               # Returns only barcode e.g. ACTGACTG
                read_barcode_distance = {barcode: Levenshtein.distance(read_barcode,barcode) for barcode in exp.picked_barcodes}
                if sum([x <= int(args.mismatch) for x in read_barcode_distance.values()]) == 0:
                    # Spacer hit but no barcode match
                    statistics["no_barcode_match"] += 1
                    continue
                if sum([x <= args.mismatch for x in read_barcode_distance.values()]) > 1:
                    # Spacer hit but multiple barcode matches
                    statistics["multiple_barcode_matches"] += 1
                    continue

                hit_barcode = min(read_barcode_distance,key=read_barcode_distance.get)

                if exp.single_cell:
                    read2 = extract_cell_barcode(read2, spacer_hit + len(args.pattern))     # The cell barcode is 16bp long and is positioned after the spacer
            
            elif not spacer_hit and not MeA_hit:
                statistics["no_spacer_found"] += 1
                # No hit, no spacer not nothing found
                if args.report_no_hit:
                    hit_barcode = 'no_spacer'
                else:
                    continue        
            
            # Now continue in the loop
            if len(read2.sequence) < 16:
                    statistics["too_short_read"] += 1
                    continue
            
            statistics[hit_barcode] += 1
            if hit_barcode in exp.picked_barcodes:
                # Write the outputs
                exp.out_stack[hit_barcode]['R1'].write('{}\n'.format(str(read1)))
                exp.out_stack[hit_barcode]['R3'].write('{}\n'.format(str(read3)))
                if args.single_cell:
                    exp.out_stack[hit_barcode]['R2'].write('{}\n'.format(str(read2)))

                
                           


    # Write the statistics file
    with open("{0}/{1}_statistics.yaml".format(exp.out_prefix,exp.name), 'w') as f:
        yaml.dump(statistics, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="DESCRIPTION: \n\nThis script demultiplexes Nano-CT sequencing data by extracting and matching modality barcodes from the R2 read based on a specified spacer sequence.\nThe script supports both bulk and single-cell data and writes sorted reads into separate output files for each detected barcode""" + 
                                     """The script works with a standard 36-8-48-36 read structure but may also be compatible with other setups where the spacer sequence is similarly arranged.\n""" + 
                                     """AATGATACGGCGACCACCGAGATCTACAC-NNNNNNNNNNNNNNNN-TCGTCGGCAGCGTCTCCACGC-NNNNNNNN-GCGATCGAGGACGGCAGATGTGTATAAGAGACAG"""+
                                     """            P5                |  sc-barcode   |  Linker sequence   | Modality |        Mosaic end               \n """ + 
                                     """""" +
                                     """ Note: If demultiplexing multiple lanes, run for each lane separetely and then merge the output files before or after alignment""",
                                     usage=""
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --single_cell --barcode ATAGAGGC                      # One specific barcode from single-cell data " 
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --single_cell --barcode ATAGAGGC TATAGCCT             # Two specific barcodes from single-cell data "
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --single_cell --Nbarcodes 3                           # Top 3 barcodes from single-cell data without specifying the barcodes - use carefuly and double check"
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --Nbarcodes 3                                         # Top 3 barcodes from bulk data ", 
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input',
                        required=True,
                        type=str,
                        nargs='+',
                        help='path to input R1,R2,R3 .fastq.gz files [3 files required]')

    parser.add_argument('-o', '--out_prefix',
                        type=str,
                        required=True,
                        help='Prefix to where to put the output files; Diretory will be created')

    parser.add_argument('-p', '--pattern',
                        type=str,
                        default="GCGTGGAGACGCTGCCGACGA",
                        help='Pattern that follows the antibody barcode \n \
                                  (Default: %(default)s)')

    parser.add_argument('--single_cell',
                        default=False,
                        action='store_true',
                        help='Data is single cell CUT&Tag (Default: %(default)s)')

    parser.add_argument('--name',
                        type=str,
                        default=None,
                        help='Custom name for the experiment (Default: Autodetect from filename)')

    parser.add_argument('--mismatch',
                        type=int,
                        default=1,
                        help='Maximum mismatches for sample barcode (Default: %(default)s)')

    parser.add_argument('--Nbarcodes',
                        type=int,
                        default=10,
                        help='Number of barcodes in experiment (Default: %(default)s)')

    parser.add_argument('--barcode',
                        type=str,
                        nargs="+",
                        default='None',
                        help='Specific barcode to be extracted [e.g. ATAGAGGC] (Default: All barcodes [see --Nbarcodes])')

    parser.add_argument('--no_barcode_seq', type=str, 
                        default='GTGTAGATCTCGGTGGTCGCCGTATCATT', 
                        help='Sequence indicating unbarcoded reads')
    
    parser.add_argument('--report_MeA', 
                        action='store_true', 
                        help='Include reads with no barcode and standard MeA sequence in the output')
    
    parser.add_argument('--report_no_hit', 
                        action='store_true', 
                        help='Include reads with no barcode or spacer whatsoever hit in the output')




    args = parser.parse_args()
    log("Starting debarcode.py script ")
    log("Input files: \n{}".format("".join(["    " + i + "\n" for i in args.input])))
    if args.barcode != "None":
        log("Provided barcodes to demultiplex: \n{}".format(args.barcode))
    log("Output prefix: {}/".format(args.out_prefix))
    
    main(args)
