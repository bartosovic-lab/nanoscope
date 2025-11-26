import argparse
from asyncio import subprocess
import yaml
import glob
import gzip
import subprocess
import re


parser = argparse.ArgumentParser(description='Summarise Nanoscope pipeline results.')
parser.add_argument('--input_config', required=True, help='Path to the input configuration file.')
# parser.add_argument('--output_summary', required=True, help='Path to the output summary file.')
args = parser.parse_args()


def main(args):
    with open(args.input_config, 'r') as f:
        config = yaml.safe_load(f)
    
    for barcode in config['barcodes']:
        modality  = config['barcodes'][barcode]
        fastq_afeter_demux_folder = f"{config['name']}/fastq_debarcoded/barcode_{barcode}/"
        all_fastq_files  = glob.glob(fastq_afeter_demux_folder + "/*R1*.fastq.gz")
        for file in all_fastq_files:
            with gzip.open(file, 'rt') as fq:
                read_count = sum(1 for line in fq) // 4
                print (f"{file}\t{barcode}\t{modality}\tfastq\tdemultiplexed\t{read_count}")
        bam_file       = f"{config['name']}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam"
        fragments_file = f"{config['name']}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz"
        cmd = ["samtools", "view", "-@", "8","-f","66", "-F", "2308", "-q", "30", "-c", bam_file]
        bam_nreads = int(subprocess.check_output(cmd, text=True).strip())
        # print(bam_nreads)
        print(f'{bam_file}\t{barcode}\t{modality}\tbam\tunique\t{bam_nreads}')

        with gzip.open(fragments_file, 'rt') as f:
            fragment_count = sum(1 for line in f if not line.startswith('#'))
        print(f'{fragments_file}\t{barcode}\t{modality}\tfragments\tcount\t{fragment_count}')

if __name__ == "__main__":
    main(args)