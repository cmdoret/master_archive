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
                                 position in another contig of a different \
                                 assembly. This is done via exact matching of \
                                 the region around the query position. This \
                                 program will send the original contig name and\
                                  position to stdout.')
parser.add_argument('--pos1', type=str, nargs='?', default=sys.stdin,
                    help='Chromosome name and basepair position in the ordered \
                    assembly. Should be in the form: Chr,bp where Chr and bp \
                    are a string and integer, respectively. Read from stdin \
                    by default.')
parser.add_argument('ref1', type=str, help='Path to the fasta file of the \
                    first assembly, where \
                    the position is known.')
parser.add_argument('ref2', type=str, help='Path to the fasta file of the \
                    second assembly, where the position is unknown.')
parser.add_argument('--region_size', type=int, default=1000, help='Size of the \
                    region to look up for exact matching in the original \
                    assembly. Default: 1000bp')

args = parser.parse_args()


def eprint(*args, **kwargs):
    # This function prints its arguments to the standard error
    print(*args, file=sys.stderr, **kwargs)

# Handling query position, whether it is from CL arg or stdin
try:
    pos1 = args.pos1.read()
except AttributeError:
    pos1 = args.pos1

# Splitting query into chromosome and basepair
if len(pos1.split(',')):
   Chr = pos1.split(',')[0]
   try:
       bp = int(pos1.split(',')[1])
   except ValueError:
       eprint("Error: bp must be an integer number representing the \
              genomic position in basepairs within the chromosome.")
else:
    eprint("Error: Query position must be in the form 'Chr,bp'.")

for chrom in SeqIO.parse(args.ref1, "fasta"):
    # Iterating over chromosomes in the ordered assembly
    if chrom.id == Chr:
        # When in correct chromosome, subset a region of defined size centered
        # around queried position
        shift_sub = [0,len(chrom.seq)]
        start_sub = bp - args.region_size / 2
        end_sub = bp + args.region_size / 2
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
                offset = bp
            if shift_sub[1] < 0:
                # If end was larger than the length of chromosome
                end_sub = len(chrom.seq)  # Shifting end to length of chrom.
                start_sub += shift_sub[1]  # Shifting start to the left
                # Offset from start to query position
                offset = region_size / 2 - shift_sub [1]
        else:
            # Region does not fit in chromosome, trim it
            offset = bp
            start_sub = 0; end_sub = len(chrom.seq)
            eprint("Warning: The lookup region did not fit inside the \
chromosome and has been trimmed down to {0} bp.".format(len(chrom.seq)))
        lookup_seq = chrom.seq[start_sub:end_sub]
        break

seq_match=[]
for contig in SeqIO.parse(args.ref2, "fasta"):
    # Iterate over all contigs in original assembly
    if contig.seq.find(lookup_seq) >= 0:
        # if the subsetted region is found in the contig, record contig name
        # and bp position corresponding to the center of the region
        seq_match.append((str(contig.id), str(contig.seq.find(lookup_seq) + offset)))

if len(seq_match) > 1:
    for match in seq_match:
        print(' '.join(match))
    eprint("Warning: {0} occurences of the lookup region were found in the\
 original assembly, you should increase --region_size.".format(len(seq_match)))
elif len(seq_match) == 1:
    print(','.join(seq_match[0]))
else:
    eprint("Error: Query position not found in reference.")
