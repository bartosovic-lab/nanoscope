import gzip
import tempfile
from os import system

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def open_fastq(file):
    return gzip.open(file, "rt") if file.endswith(".gz") else open(file, "r")

def concat_reads(r1_file, r2_file, output_file):
    handle_r1 = open_fastq(r1_file)
    handle_r2 = open_fastq(r2_file)
    output_handle = gzip.open(output_file, "wt") if output_file.endswith(".gz") else open(output_file, "w")

    r1_iter = SeqIO.parse(handle_r1, "fastq")
    r2_iter = SeqIO.parse(handle_r2, "fastq")

    for rec1, rec2 in zip(r1_iter, r2_iter):
        assert rec1.id.split()[0] == rec2.id.split()[0], f"ID mismatch: {rec1.id} != {rec2.id}"

        combined_seq = rec1.seq + rec2.seq
        combined_qual = rec1.letter_annotations["phred_quality"] + rec2.letter_annotations["phred_quality"]

        combined_record = SeqRecord(
            combined_seq,
            id=rec1.id,
            description="",
            letter_annotations={"phred_quality": combined_qual}
        )

        SeqIO.write(combined_record, output_handle, "fastq")

    handle_r1.close()
    handle_r2.close()
    output_handle.close()

def run_je_demultiplex(combined_file, R2_file, barcodes_table,output_folder="demux_out"):
    command = [
        "je debarcode",
        "F={}".format(combined_file),
        "F={}".format(R2_file),
        "BF={}".format(barcodes_table),
        "RL='<SAMPLE1:x>'",
        "RL='<BARCODE1:8><SAMPLE3:x>'",
        "O={}".format(output_folder),
        "GZ=true",
        "UN=true",
        "ADD=true",
        "MM=1",
        "MMD=1",
        "CLIP=true"
    ]
    print("Running command:", " ".join(command))
    system(" ".join(command) if isinstance(command, list) else command)

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Concatenate paired-end reads from two FASTQ files.")
    parser.add_argument("-1", "--r1_file", help="Input FASTQ file for read 1 (can be gzipped)")
    parser.add_argument("-2", "--r2_file", help="Input FASTQ file for read 2 (can be gzipped)")
    parser.add_argument("-3", "--r3_file", help="Input FASTQ file for read 3 (can be gzipped)")
    parser.add_argument("-b", "--barcodes_table", help="Barcodes table for demultiplexing")
    parser.add_argument("-o", "--output_file", help="Output FASTQ file for concatenated reads")
    args = parser.parse_args()

    with tempfile.NamedTemporaryFile(delete=True, suffix=".fastq.gz") as temp_file:
        concat_reads(args.r1_file, args.r3_file, temp_file.name)
        run_je_demultiplex(temp_file.name, args.r2_file, barcodes_table=args.barcodes_table)

if __name__ == "__main__":
    main()