# This is a simple python script to obtain the position and contig in a
# reference fasta file, given the position and contig in another fasta file.
# Both fasta files must be provided, and exact match search is performed in a
# region of given length.
# Cyril Matthey-Doret
# 19.10.2017

from __future__ import print_function  # Use python3 print function
import sys  # Will allow to print messages to stderr
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description = 'This script allows to retrieve \
                                 the position in an original contig from a \
                                 position in a chromosome of the ordered \
                                 assembly obtained with the linkage map. This \
                                 program will send the original contig name and\
                                  position to stdout.')
parser.add_argument('Chr', type=str,help='Chromosome name in the ordered \
                    assembly.')
parser.add_argument('bp', type=int, help='Integer position, in basepair, \
                    within the chromosome of the ordered assembly.')
parser.add_argument('ord_ref', type=str, help='Path to the ordered assembly')
parser.add_argument('orig_ref', type=str, help='Path to the original assembly.')
parser.add_argument('--region_size', type=int, default=1000, help='Size of the \
                    region to look up for exact matching in the original \
                    assembly. Default: 1000bp')


args = parser.parse_args()

def eprint(*args, **kwargs):
    # This function prints its arguments to the standard error
    print(*args, file=sys.stderr, **kwargs)

## TEST

for chrom in SeqIO.parse(args.ord_ref, "fasta"):
    # Iterating over chromosomes in the ordered assembly
    if chrom.id == args.Chr:
        # When in correct chromosome, subset a region of defined size centered
        # around queried position
        shift_sub = [0,len(chrom.seq)]
        start_sub = args.bp - args.region_size / 2
        end_sub = args.bp + args.region_size / 2
        offset = args.region_size / 2
        if args.region_size <= len(chrom.seq):
            # Region fits in chromosome, shift if needed
            shift_sub[0] -= start_sub
            shift_sub[1] -= end_sub
            if shift_sub[0] > 0:
                # If start had a negative value
                start_sub = 0  # Shifting start to zero
                end_sub += shift_sub[0]  # Shifting  end to the right accordingly
                # Offset from start to query position
                offset = args.bp
            if shift_sub[1] < 0:
                # If end was larger than the length of chromosome
                end_sub = len(chrom.seq)  # Shifting end to length of chrom.
                start_sub += shift_sub[1]  # Shifting start to the left
                # Offset from start to query position
                offset = region_size / 2 - shift_sub [1]
        else:
            # Region does not fit in chromosome, trim it
            offset = args.bp
            start_sub = 0; end_sub = len(chrom.seq)
            eprint("Warning: The lookup region did not fit inside the \
chromosome and has been trimmed down to {0} bp.".format(len(chrom.seq)))
        lookup_seq = chrom.seq[start_sub:end_sub]
        break

unord_match=[]
for contig in SeqIO.parse(args.orig_ref, "fasta"):
    # Iterate over all contigs in original assembly
    if contig.seq.find(lookup_seq) >= 0:
        # if the subsetted region is found in the contig, record contig name
        # and bp position corresponding to the center of the region
        unord_match.append((str(contig.id), str(contig.seq.find(lookup_seq) + offset)))

if len(unord_match) > 1:
    for match in unord_match:
        print(' '.join(match))
    eprint("Warning: {0} occurences of the lookup region were found in the\
 original assembly, you should increase --region_size.".format(len(unord_match)))
else:
    print(' '.join(unord_match[0]))
