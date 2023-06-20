import pysam
import sys
import json

# Marek Bartosovic 13/06/2023
# Script which removes linear amplification duplicates from possorted_bam.bam file from cellranger
# Script keeps and reports PCR duplicates as these are handled fine in downstream analysis

# Warning: Input possorted_bam.bam file needs to be re-sorted by read name otherwise the script won't work
# This can be done using following command:
# samtools sort -@ 10 -n -o namesorted_bam.bam possorted_bam.bam


# Input
samfile = pysam.AlignmentFile(sys.argv[1], "rb")

#Output
sam_out = pysam.AlignmentFile(sys.argv[2], "wb", header=samfile.header)

# Stats output:
stats_out = sys.argv[2] + "_stats.txt"

reads_unique = {}
stats = {
    'not mapped': 0,
    'unique': 0,
    'not proper': 0,
    'LA duplicates': 0,
    'PCR duplicates': 0,
}

n = 0
# n_target = 1000000  # test for 1 milion reads
# n_target = 10000000 # test for 10 milion reads

read1 = False
read2 = False

for line in samfile:
    n += 1
    if n % 100000 == 0:
        sys.stderr.write('*** {} lines processed\n'.format(n))

    # # Only take the first n_target reads # Debugging
    # if n == n_target:
    #     # print(reads_unique)
    #     print(stats)
    #     break

    # File needs to be namesorted - read a read pair
    read1 = line            # R1
    read2 = next(samfile)   # R2

    # Some assertion to make sure we are in namesorted file and read1 and read2 follow each other
    assert read1.query_name == read2.query_name, "Reads name do not match for R1 and R2 in: \n" + str(read1) + "\n" + str(read2) + '\n *** Perhaps the bam file is not namesorted ? ***'
    assert read1.is_read1
    assert not read2.is_read1

    if read1.is_unmapped or read2.is_unmapped:
        stats['not mapped'] += 1
        continue

    if not read1.is_proper_pair or not read2.is_proper_pair:
        stats['not proper'] += 1
        continue

    cell_barcode = read1.get_tag('CR')
    read1_position = '{}_{}_{}'.format(read1.reference_name, + read1.reference_start, cell_barcode)
    read2_position = '{}_{}_{}'.format(read1.next_reference_name, + read1.next_reference_start, cell_barcode)

    try:  # Creates new fw read entry in the dictionary if not there. # Dict is a lot faster if indexed only by string and not whole AlignedSegment object
        reads_unique[read1_position]
    except KeyError:
        reads_unique[read1_position] = []

    # Case1 read2_position list is empty == new read is observed             == new unique read
    if reads_unique[read1_position] == []:
        reads_unique[read1_position].append(read2_position)
        stats['unique'] += 1
        sam_out.write(read1)
        sam_out.write(read2)
        continue

    # Case read2_position list is not empty and the read2_position is not in list    == LA duplicate
    if read2_position not in reads_unique[read1_position]:
        reads_unique[read1_position].append(read2_position)
        stats['LA duplicates'] += 1
        # LA duplicates are not reported in the final bam file
        continue

    # Case read2_position list is not empty and read2_position is in the list        == PCR duplicate
    if read2_position in reads_unique[read1_position]:
        stats['PCR duplicates'] += 1
        # PCR duplicates can be reported as they are handled in downstream analysis
        sam_out.write(read1)
        sam_out.write(read2)
        continue

    # Sanity check if something escapes
    sys.exit('found exception')

with open(stats_out,'w') as f:
    json.dump(stats,f)
