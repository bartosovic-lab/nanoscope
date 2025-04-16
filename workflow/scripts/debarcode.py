import os
import sys
import argparse
import gzip
import time
import yaml
import regex
import Levenshtein
from glob import glob
from contextlib import ExitStack
from collections import defaultdict
from re import split
from pysam import FastxFile


def log(msg):
    sys.stderr.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}: {msg}\n")

class BcdCT:
    def __init__(self, args):
        self.args = args
        self._detect_input(args.input)      # <- moved up here!
        self._detect_reads()

        self.single_cell = args.single_cell
        self.out_reads = ['R1', 'R2', 'R3'] if self.single_cell else ['R1', 'R3']
        self.name = args.name or self._infer_name()

        self.provided_barcodes = args.barcode if args.barcode != "None" else None
        if self.provided_barcodes:
            self._validate_provided_barcodes()

        self._autodetect_barcodes()
        self._filter_valid_barcodes()
        self._prepare_output_paths()

        print(f"Provided barcodes: {self.provided_barcodes}")
        print(f"Top barcodes (n={args.Nbarcodes}): {self.top_barcodes}")
        print(f"Picked barcodes: {self.picked_barcodes}")

    def _detect_input(self, paths):
        paths = [os.path.abspath(p) for p in paths]
        if len(paths) == 1 and os.path.isdir(paths[0]):
            self.input_files = glob(os.path.join(paths[0], '*.f*q.gz'))
        elif all(p.endswith(('.fastq.gz', '.fq.gz')) for p in paths):
            self.input_files = paths
        else:
            sys.exit("*** Error: Input must be a folder or .fastq.gz/.fq.gz files with _R1_, _R2_, _R3_ ***")

    def _detect_reads(self):
        self.path_in = {r: [f for f in self.input_files if f"_{r}_" in f] for r in ['R1', 'R2', 'R3']}
        missing = [r for r, f in self.path_in.items() if len(f) != 1]
        if missing:
            sys.exit(f"*** Error: Exactly one file for each of R1, R2, R3 required. Missing: {missing} ***")
        self.path_in = {r: f[0] for r, f in self.path_in.items()}

    def _infer_name(self):
        prefixes = {split("_R[0-9]_", f)[0] for f in self.input_files}
        if len(prefixes) > 1:
            sys.exit("*** Error: Input files do not share a common prefix. Use --name. ***")
        return os.path.basename(prefixes.pop())

    def _validate_provided_barcodes(self):
        allowed_special = {"MeA"}
        valid = []
        for b in self.provided_barcodes:
            if all(c in "ATCGN" for c in b) or b in allowed_special:
                valid.append(b)
            else:
                sys.stderr.write(f"*** Warning: Invalid barcode '{b}' ignored. Must be A/T/C/G/N or allowed keyword like 'MeA'. ***\n")
        self.provided_barcodes = valid

    def _filter_valid_barcodes(self):
        allowed_special = {"MeA", "no_hit"}
        self.valid_barcodes = {
            b for b in self.picked_barcodes
            if all(c in "ATCGN" for c in b) or b in allowed_special
        }

        if not self.valid_barcodes:
            sys.exit("*** Error: No valid barcodes found. ***")

    def _autodetect_barcodes(self):
        log("Detecting barcodes...")
        counts = defaultdict(int)
        for i, (_, r2, _) in enumerate(self._read_triplets()):
            hit = find_seq(self.args.pattern, r2.sequence, 0)
            if hit is not None:                
                barcode = revcompl(r2.sequence[hit - 8:hit])
                counts[barcode] += 1
            elif find_seq(self.args.no_barcode_seq, r2.sequence, 0):
                counts['MeA'] += 1
            else:
                counts['no_hit'] += 1
            if i >= 50000:
                break
        self.top_barcodes = dict(sorted(counts.items(), key=lambda x: x[1], reverse=True)[:self.args.Nbarcodes])
        self.picked_barcodes = (
            list(self.provided_barcodes) if self.provided_barcodes else list(self.top_barcodes.keys())
        )
        if self.args.report_no_hit and 'no_hit' not in self.picked_barcodes:
            self.picked_barcodes.append('no_hit')
        if self.args.report_MeA and 'MeA' not in self.picked_barcodes:
            self.picked_barcodes.append('MeA')

    

    def _prepare_output_paths(self):
        log("Preparing output paths...")
        self.path_out = {
            b: {
                r: os.path.join(self.args.out_prefix, f"barcode_{b}", os.path.basename(self.path_in[r]).replace(split('_S[0-9]+_', os.path.basename(self.path_in[r]))[0].strip("_"), self.name))
                for r in self.out_reads
            } for b in self.picked_barcodes
        }

    def _read_triplets(self):
        with FastxFile(self.path_in['R1']) as r1, FastxFile(self.path_in['R2']) as r2, FastxFile(self.path_in['R3']) as r3:
            yield from zip(r1, r2, r3)

    def create_output_handles(self, stack):
        for b in self.picked_barcodes:
            os.makedirs(os.path.join(self.args.out_prefix, f"barcode_{b}"), exist_ok=True)
        self.out_stack = {
            b: {r: stack.enter_context(gzip.open(self.path_out[b][r], 'wt')) for r in self.out_reads}
            for b in self.picked_barcodes
        }


def revcompl(seq):
    return seq.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def find_seq(pattern, seq, nmismatch=2):
    for n in range(nmismatch + 1):
        matches = list(regex.finditer(f'({pattern}){{e<={n}}}', seq))
        if len(matches) == 1:
            return matches[0].start()
    return None


def extract_cell_barcode(read, idx, length=16):
    read.sequence = read.sequence[idx:idx + length]
    read.quality = read.quality[idx:idx + length]
    return read


def main(args):
    exp = BcdCT(args)
    stats = defaultdict(int)

    with ExitStack() as stack:
        exp.create_output_handles(stack)
        for i, (r1, r2, r3) in enumerate(exp._read_triplets()):
            if i % 5_000_000 == 0:
                log(f"{i:,} reads processed")

            if r1.name != r2.name or r2.name != r3.name:
                continue

            hit = find_seq(args.pattern, r2.sequence, 2)
            matched_barcode = None

            if hit is None:
                if find_seq(args.no_barcode_seq, r2.sequence, 0):
                    stats["MeA"] += 1
                    matched_barcode = 'MeA'
                    hit = 0
                else:
                    stats["no_hit"] += 1
                    matched_barcode = 'no_hit'
                    continue

            if hit:
                barcode = revcompl(r2.sequence[hit - 8:hit])
                close = {b: Levenshtein.distance(barcode, b) for b in exp.valid_barcodes if Levenshtein.distance(barcode, b) <= args.mismatch}
                if not close:
                    stats["no_barcode_match"] += 1
                    continue
                if len(close) > 1:
                    stats["ambiguous_barcode_matches"] += 1
                matched_barcode = min(close, key=close.get)

            if matched_barcode in exp.picked_barcodes:
                if exp.single_cell:
                    if matched_barcode == "MeA":
                        if len(r2.sequence) >= 16:
                            r2 = extract_cell_barcode(r2, 0)
                        else:
                            stats["too_short_read"] += 1
                            continue
                    elif hit is not None:
                        start = hit + len(exp.args.pattern)
                        if len(r2.sequence) >= start + 16:
                            r2 = extract_cell_barcode(r2, start)
                        else:
                            stats["too_short_read"] += 1
                            continue
                    else:
                        stats["too_short_read"] += 1
                        continue
                    r2.sequence = revcompl(r2.sequence)
                    exp.out_stack[matched_barcode]['R2'].write(f"{r2}\n")
                exp.out_stack[matched_barcode]['R1'].write(f"{r1}\n")
                exp.out_stack[matched_barcode]['R3'].write(f"{r3}\n")

    with open(os.path.join(args.out_prefix, f"{exp.name}_statistics.yaml"), 'w') as f:
        yaml.dump(dict(stats), f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Demultiplex Nano-CT sequencing data using modality barcodes in R2 read.",
                                    formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', type=str, nargs='+', required=True, help='Input R1, R2, R3 .fastq.gz files or input folder')
    parser.add_argument('-o', '--out_prefix', type=str, required=True, help='Output folder prefix')
    parser.add_argument('-p', '--pattern', type=str, default="GCGTGGAGACGCTGCCGACGA", help='Pattern following the barcode')
    parser.add_argument('--single_cell', action='store_true', help='Enable single-cell mode')
    parser.add_argument('--name', type=str, help='Custom experiment name (default: inferred)')
    parser.add_argument('--mismatch', type=int, default=2, help='Max barcode mismatches')
    parser.add_argument('--Nbarcodes', type=int, default=6, help='Number of top barcodes to select')
    parser.add_argument('--barcode', nargs='+', default='None', help='Specific barcodes to extract')
    parser.add_argument('--no_barcode_seq', type=str, default='GTGTAGATCTCGGTGGTCGCCGTATCATTAAA', help='Sequence indicating unbarcoded reads')
    parser.add_argument('--report_MeA', action='store_true', help='Include short MeA reads in output')
    parser.add_argument('--report_no_hit', action='store_true', help='Include reads with no barcode hit')

    main(parser.parse_args())