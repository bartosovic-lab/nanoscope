from Bio import SeqIO
import argparse
import gzip
import pysam
import os


parser = argparse.ArgumentParser(description="Process R2 barcodes coming from multiplexed, multimodal Tn5 experiments")
parser.add_argument("--input_fastq", required=True, help="Input FASTQ file of format" + 
                                    """AATGATACGGCGACCACCGAGATCTACAC-NNNNNNNNNNNNNNNN-TCGTCGGCAGCGTCTCCACGC-NNNNNNNN-GCGATCGAGGACGGCAGATGTGTATAAGAGACAG""" +
                                    """            P5                |  sc-barcode   |  Linker sequence   | Modality |        Mosaic end            \n """ + 
                                    """                                                                               <- sequencing primer             """)
parser.add_argument("--output_modality_suffix", default="_modality", help="Suffix for output modality FASTQ files")
parser.add_argument("--output_scbarcode_suffix", default="_singlecell", help="Suffix for output sc-barcode FASTQ files")
parser.add_argument("--output_folder", default=".", help="Output folder for processed FASTQ files")
args = parser.parse_args()

os.makedirs(args.output_folder, exist_ok=True)

def main():
    outfile_modality = '{folder}/{fastq_in}{suffix}'.format(folder=args.output_folder, \
                                                            fastq_in=os.path.basename(args.input_fastq).replace('.fastq.gz',''), \
                                                            suffix=args.output_modality_suffix + '.fastq.gz')
    outfile_scbarcode = '{folder}/{fastq_in}{suffix}'.format(folder=args.output_folder, \
                                                            fastq_in=os.path.basename(args.input_fastq).replace('.fastq.gz',''), \
                                                            suffix=args.output_scbarcode_suffix + '.fastq.gz')

    print(args.input_fastq)
    print(outfile_modality,outfile_scbarcode)
    with gzip.open(outfile_modality, "wt") as out_mod_fh, gzip.open(outfile_scbarcode, "wt") as out_sc_fh:
        for read in pysam.FastxFile(args.input_fastq):
            # Create new read objects for modality and sc-barcode
            modality = pysam.FastxRecord(
                name=read.name,
                sequence=read.sequence[:8],  # Extract modality (8 bases after linker)
                quality=read.quality[:8])
            sc_barcode = pysam.FastxRecord(
                name=read.name,
                sequence=read.sequence[8+21:8+21+16],  # Extract sc-barcode (
                quality=read.quality[8+21:8+21+16])
            # print(modality)
            # print(out_sc_barcode)
            out_mod_fh.write(str(modality) + "\n")
            out_sc_fh.write(str(sc_barcode) + "\n")

main()
