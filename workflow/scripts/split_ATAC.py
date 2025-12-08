import argparse
import pysam
import regex
import os
import gzip
import copy


parser = argparse.ArgumentParser(description="Simple debarcoding script")
parser.add_argument("-i","--input", help="Input FASTQ file R1 R2 R3", nargs= '+', required=True)
parser.add_argument("--ATAC", default="ATAC", help="ATAC FASTQ output directory")
parser.add_argument("--nonATAC", default="nonATAC", help="non-ATAC FASTQ output directory")
parser.add_argument("--noMatch", default = "noMatch", help="no-match FASTQ output directory")
args = parser.parse_args()

# Pattern with up to 2 mismatches allowed
ATAC_pattern   = regex.compile(r'(GTGTAGATCTCGGTGGTCGCCGTATCATTAAA){e<=2}', regex.IGNORECASE)
nanoCT_pattern = regex.compile(r'(GCGTGGAGACGCTGCCGACGA){e<=2}', regex.IGNORECASE)

n = 0
n_final = 10000  # for testing

def parse_input(input_files):
    types = ['R1', 'R2', 'R3']
    input_dict = {}
    for t in types:
        input_dict[t] = [
            f for f in input_files
            if t in os.path.basename(f)
            and (f.endswith('.fastq') or f.endswith('.fq') or
                 f.endswith('.fastq.gz') or f.endswith('.fq.gz'))
        ]
        if len(input_dict[t]) > 1:
            raise ValueError(f"Multiple {t} files provided: {input_dict[t]}")
        elif len(input_dict[t]) == 0:
            raise ValueError(f"No {t} file provided.")
        else:
            input_dict[t] = input_dict[t][0]
    return input_dict

def smart_open(path, mode="rt"):
    """Open plain or gzipped text file."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    else:
        return open(path, mode)

# Create ATAC/nonATAC directories if they don't exist
os.makedirs(args.ATAC, exist_ok=True)
os.makedirs(args.nonATAC, exist_ok=True)
os.makedirs(args.noMatch, exist_ok=True)
input_dict = parse_input(args.input)

with pysam.FastxFile(input_dict['R1']) as R1_file, \
     pysam.FastxFile(input_dict['R2']) as R2_file, \
     pysam.FastxFile(input_dict['R3']) as R3_file, \
     smart_open(f"{args.ATAC}/{os.path.basename(input_dict['R1'])}", "wt") as ATAC_R1_file, \
     smart_open(f"{args.ATAC}/{os.path.basename(input_dict['R2'])}", "wt") as ATAC_R2_file, \
     smart_open(f"{args.ATAC}/{os.path.basename(input_dict['R3'])}", "wt") as ATAC_R3_file, \
     smart_open(f"{args.nonATAC}/{os.path.basename(input_dict['R1'])}", "wt") as nonATAC_R1_file, \
     smart_open(f"{args.nonATAC}/{os.path.basename(input_dict['R2'])}", "wt") as nonATAC_R2_file, \
     smart_open(f"{args.nonATAC}/{os.path.basename(input_dict['R3'])}", "wt") as nonATAC_R3_file, \
     smart_open(f"{args.nonATAC}/{os.path.basename(input_dict['R2'])}_modality.gz", "wt") as nonATAC_R2_modality_file, \
     smart_open(f"{args.noMatch}/{os.path.basename(input_dict['R1'])}", "wt") as noMatch_R1_file, \
     smart_open(f"{args.noMatch}/{os.path.basename(input_dict['R2'])}", "wt") as noMatch_R2_file, \
     smart_open(f"{args.noMatch}/{os.path.basename(input_dict['R3'])}", "wt") as noMatch_R3_file:
    for R1_read, R2_read, R3_read in zip(R1_file, R2_file, R3_file):
        n += 1
        if n >= n_final:
            break
        match_ATAC   = ATAC_pattern.search(R2_read.sequence)
        match_nanoCT = nanoCT_pattern.search(R2_read.sequence)
        if (match_nanoCT and match_ATAC) or (not match_nanoCT and not match_ATAC):
            # Write to no-match files
            noMatch_R1_file.write(str(R1_read) + "\n")
            noMatch_R2_file.write(str(R2_read) + "\n")
            noMatch_R3_file.write(str(R3_read) + "\n")
            continue
        if match_ATAC:
            # Trim R2 read to first 16 bases (cell barcode)
            R2_read.sequence = R2_read.sequence[:16]
            R2_read.quality  = R2_read.quality[:16]
            # Write to ATAC files
            ATAC_R1_file.write(str(R1_read) + "\n")
            ATAC_R2_file.write(str(R2_read) + "\n")
            ATAC_R3_file.write(str(R3_read) + "\n")
            continue
        if match_nanoCT:
            R2_read_singlecell = copy.copy(R2_read)
            R2_read_modality   = copy.copy(R2_read)
            
            # Modality is bases 1-8
            R2_read_modality.sequence = R2_read_modality.sequence[:8]
            R2_read_modality.quality  = R2_read_modality.quality[:8]

            # GCGTGGAGACGCTGCCGACGA follows and then 16 bp single-cell index
            R2_read_singlecell.sequence = R2_read_singlecell.sequence[29:45]
            R2_read_singlecell.quality  = R2_read_singlecell.quality[29:45]

            # Write to non-ATAC files
            nonATAC_R1_file.write(str(R1_read) + "\n")
            nonATAC_R2_file.write(str(R2_read_singlecell) + "\n")
            nonATAC_R2_modality_file.write(str(R2_read_modality) + "\n")
            nonATAC_R3_file.write(str(R3_read) + "\n")
            continue
        else:
            raise ValueError("This should not happen")