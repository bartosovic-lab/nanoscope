import sys
import argparse
import gzip

parser = argparse.ArgumentParser(description='Filter cellranger gtf file and return genebody and promter gtf file and separate file with gene names')
parser.add_argument('-i','--input', type=str, help='GTF file')
parser.add_argument('-o','--output', type=str, help='Output GTF file')
parser.add_argument('-n','--gene_names', type=str, help='Output gene names file')
args = parser.parse_args()

def parse_9th_column(col9):
    col9 = col9.split(";")
    col9 = {x.split()[0]: x.split()[1].replace('"', '') for x in col9 if x}
    return col9

def extend_promoter(line,by=2000):
    if line[6] == "+":
        line[3] = int(line[3]) - by
        # Prevent negative values
        if line[3] < 0:
            line[3] = 1
    else:
        line[4] = int(line[4]) + by
    line[3] = str(line[3])
    line[4] = str(line[4])
    return line

def check_gene_name(col9):
    if "gene_name" in col9:
        gname_col = "gene_name"
    elif "gene_id" in col9:
        gname_col = "gene_id"
    else:
        gname_col = col9.keys()[0]
    return gname_col


def main(args):
    if args.input.endswith(".gz"):
        gtf = gzip.open(args.input, "rt")
    else:
        gtf = open(args.input, "r")
    with open(args.output, "w") as out, open(args.gene_names, "w") as gene_names:
        for line in gtf:
            if line.startswith("#"):
                out.write(line)
                continue

            line = line.strip().split("\t")
            if line[2] == "gene":
                col9 = parse_9th_column(line[8])
                line = extend_promoter(line,by=2000)
                out.write("\t".join(line) + "\n")
                try:    
                    gene_names.write(col9["gene_name"] + "\n")
                except:
                    gene_names.write(col9["gene_id"] + "\n")

if __name__ == "__main__":
    main(args)
    