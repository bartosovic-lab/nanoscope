import os
import sys
import argparse
import gzip
from glob import glob
from contextlib import ExitStack
from collections import defaultdict
from re import split
import regex
import yaml
from pysam import FastxFile
import Levenshtein

class BcdCT:
    def __init__(self, args):
        self.args = args
        self.detect_input(args.input)
        self.detect_reads()

        self.single_cell = args.single_cell
        self.out_reads = ['R1', 'R2', 'R3'] if self.single_cell else ['R1', 'R3']
        self.name = args.name or self.autodetect_name()

        self.autodetect_barcodes()
        self.prepare_output_filenames()

    def detect_input(self, input_paths):
        input_paths = [os.path.abspath(path) for path in input_paths]
        error_message = (
            "*** Error: Wrong input files specified. The input must be a folder or paths to fastq.gz files. ***\n"
            "Expecting _R1_, _R2_, _R3_ in filenames.\n"
        )

        if len(input_paths) == 1 and os.path.isdir(input_paths[0]):
            self.input_dir = input_paths[0]
            self.input_files = glob(os.path.join(self.input_dir, '*.fastq.gz')) + \
                               glob(os.path.join(self.input_dir, '*.fq.gz'))
        elif len(input_paths) > 1:
            self.input_files = input_paths
            self.input_dir = os.path.dirname(input_paths[0])
            if not all(f.endswith(('.fastq.gz', '.fq.gz')) for f in self.input_files):
                sys.exit(error_message)
        else:
            sys.exit(error_message)

    def detect_reads(self):
        """Identify and verify R1, R2, and R3 FASTQ files from input files."""
        self.path_in = {
            read: [f for f in self.input_files if f"_{read}_" in f]
            for read in ['R1', 'R2', 'R3']
        }

        missing = [read for read in self.path_in if len(self.path_in[read]) != 1]
        if missing:
            sys.exit(f"*** Error: Must have exactly one file for each of R1, R2, and R3. Missing or ambiguous: {missing} ***\n")

        # Flatten each list to a single path
        self.path_in = {read: files[0] for read, files in self.path_in.items()}

    def prepare_output_filenames(self):
        self.path_out = {
            barcode: {
                read: os.path.join(self.args.out_prefix, f"barcode_{barcode}", os.path.basename(self.path_in[read]))
                for read in self.out_reads
            } for barcode in self.picked_barcodes
        }

        if self.args.name:
            for barcode, reads in self.path_out.items():
                for read_type, path in reads.items():
                    sample_id = split('_S[0-9]+_', os.path.basename(path))[0].strip("_")
                    self.path_out[barcode][read_type] = path.replace(sample_id, self.args.name)

    def autodetect_name(self):
        prefixes = [split("_R[0-9]_", f)[0] for f in self.path_in.values()]
        unique_prefixes = set(prefixes)
        if len(unique_prefixes) > 1:
            sys.exit("*** Error: Input files do not share the same prefix. Use --name to specify one. ***\n")
        return os.path.basename(prefixes[0])

    def in_handles(self, stack):
        return {read: stack.enter_context(FastxFile(self.path_in[read], 'r')) for read in ['R1', 'R2', 'R3']}

    def create_output_handles(self, stack):
        for barcode in self.picked_barcodes:
            os.makedirs(os.path.join(self.args.out_prefix, f"barcode_{barcode}"), exist_ok=True)
        self.out_stack = {
            barcode: {
                read: stack.enter_context(gzip.open(self.path_out[barcode][read], 'wt'))
                for read in self.out_reads
            } for barcode in self.picked_barcodes
        }

    def autodetect_barcodes(self):
        barcode_counts = defaultdict(int)
        barcode_counts['no_barcode'] = 0
        barcode_counts['no_hit'] = 0

        for n, (r1, r2, r3) in enumerate(self):
            hit = find_seq(self.args.pattern, r2.sequence, nmismatch=0)
            if not hit:
                if find_seq(self.args.no_barcode_seq, r2.sequence, nmismatch=0):
                    barcode_counts['no_barcode'] += 1
                else:
                    barcode_counts['no_hit'] += 1
            else:
                barcode = get_read_barcode(r2, hit)
                barcode_counts[barcode] += 1
            if n >= 50000:
                break

        sorted_barcodes = sorted(barcode_counts.items(), key=lambda x: x[1], reverse=True)
        top_barcodes = dict(sorted_barcodes[:self.args.Nbarcodes])

        print(f"Top barcodes (n={n}): {top_barcodes}")
        if self.args.barcode != "None":
            self.picked_barcodes = self.args.barcode
        else:
            self.picked_barcodes = list(top_barcodes.keys())
        # Always include 'no_barcode' and 'no_hit' in the outputs
        if 'no_hit' not in self.picked_barcodes:
            self.picked_barcodes.append('no_hit')
        if 'no_barcode' not in self.picked_barcodes:
            self.picked_barcodes.append('no_barcode')
        print(f"Picked barcodes: {self.picked_barcodes}")

    def __iter__(self):
        with FastxFile(self.path_in['R1']) as f1, FastxFile(self.path_in['R2']) as f2, FastxFile(self.path_in['R3']) as f3:
            yield from zip(f1, f2, f3)

# === Helper functions ===

def get_read_barcode(read, index):
    return revcompl(read.sequence[index - 8:index])

def extract_cell_barcode(read, index, pattern):
    start = index + len(pattern)
    read.sequence = read.sequence[start:start + 16]
    read.quality = read.quality[start:start + 16]
    return read

def revcompl(seq):
    rev_table = str.maketrans("ACGTN", "TGCAN")
    return seq.upper().translate(rev_table)[::-1]

def find_seq(pattern, sequence, nmismatch=2):
    for n in range(nmismatch + 1):
        matches = list(regex.finditer(f'({pattern}){{e<={n}}}', sequence))
        if len(matches) == 1:
            return matches[0].start()
        elif len(matches) > 1:
            return None
    return None

# === Main ===

def main(args):
    exp = BcdCT(args)
    stats = defaultdict(int)

    with ExitStack() as stack:
        exp.create_output_handles(stack)
        for n, (r1, r2, r3) in enumerate(exp):
            if n % 5_000_000 == 0:
                print(f"{n} reads processed", file=sys.stderr)
            if r1.name != r2.name or r2.name != r3.name:
                continue

            hit = find_seq(args.pattern, r2.sequence, nmismatch=2)
            if not hit and find_seq(args.no_barcode_seq, r2.sequence, nmismatch=0):
                hit = 0
                stats["no_barcode"] += 1

            if hit is not None:
                barcode = get_read_barcode(r2, hit)
                distances = {b: Levenshtein.distance(barcode, b) for b in exp.picked_barcodes}
                close_matches = [b for b, d in distances.items() if d <= args.mismatch]

                if not close_matches:
                    stats["no_barcode_match"] += 1
                    continue
                if len(close_matches) > 1:
                    stats["multiple_barcode_matches"] += 1
                    continue

                matched_barcode = close_matches[0]

                if exp.single_cell:
                    r2 = extract_cell_barcode(r2, hit, args.pattern)
                    if len(r2.sequence) < 16:
                        stats["too_short_read"] += 1
                        continue
                    exp.out_stack[matched_barcode]['R2'].write(f"{r2}\n")

                exp.out_stack[matched_barcode]['R1'].write(f"{r1}\n")
                exp.out_stack[matched_barcode]['R3'].write(f"{r3}\n")

    with open(os.path.join(args.out_prefix, f"{exp.name}_statistics.yaml"), 'w') as f:
        yaml.dump(dict(stats), f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="DESCRIPTION: \n\nThis script demultiplexes Nano-CT sequencing data by extracting and matching modality barcodes from the R2 read based on a specified spacer sequence.\nThe script supports both bulk and single-cell data and writes sorted reads into separate output files for each detected barcode\n""" + 
                                     """The script works with a standard 36-8-48-36 read structure but may also be compatible with other setups where the spacer sequence is similarly arranged.\n\n""" + 
                                     """AATGATACGGCGACCACCGAGATCTACAC-NNNNNNNNNNNNNNNN-TCGTCGGCAGCGTCTCCACGC-NNNNNNNN-GCGATCGAGGACGGCAGATGTGTATAAGAGACAG\n"""+
                                     """            P5                |  sc-barcode   |  Linker sequence   | Modality |        Mosaic end               \n """ + 
                                     """\n""" +
                                     """ Note: If demultiplexing multiple lanes, run for each lane separetely and then merge the output files before or after alignment\n""",
                                     usage="\n"
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --single_cell --barcode ATAGAGGC                      # One specific barcode from single-cell data \n" 
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --single_cell --barcode ATAGAGGC TATAGCCT             # Two specific barcodes from single-cell data \n"
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --single_cell --Nbarcodes 3                           # Top 3 barcodes from single-cell data without specifying the barcodes - use carefuly and double check\n"
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --Nbarcodes 3                                         # Top 3 barcodes from bulk data \n", 
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
                        default=3,
                        help='Number of barcodes in experiment (Default: %(default)s)')

    parser.add_argument('--barcode',
                        type=str,
                        nargs="+",
                        default='None',
                        help='Specific barcode to be extracted [e.g. ATAGAGGC] (Default: All barcodes [see --Nbarcodes])')
    parser.add_argument('--no_barcode_seq',default = 'GTGTAGATCTCGGTGGTCGCCGTATCATTAAA',
                        type=str,
                        help='Sequence to be used for demultiplexing of unbarcoded reads (Default: %(default)s)')
    
    args = parser.parse_args()
    print(args.barcode)
    main(args)
