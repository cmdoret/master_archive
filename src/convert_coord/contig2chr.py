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

# PARSE CL ARGS

parser = argparse.ArgumentParser(description='This script allows to retrieve \
                                 the position in an original contig from a \
                                 position in another contig of a different \
                                 assembly. This is done via exact matching of \
                                 the region around the query position. This \
                                 program sends the original contig name and\
                                  position to stdout. It will also return the \
                                  region size used and if the region has been \
                                  reversed and/or complemented.')
parser.add_argument('--pos1', type=str, nargs='?', default=sys.stdin,
                    help='Chromosome name and basepair position in the ordered \
                    assembly. Should be in the form: Chr,bp where Chr and bp \
                    are a string and integer, respectively. Read from stdin \
                    by default.')
parser.add_argument('ref1', type=str, help='Path to the fasta file of the \
                    assembly, where the position is known.')
parser.add_argument('ref2', type=str, help='Path to the fasta file of the \
                    assembly, where the position is unknown.')
parser.add_argument('--region_size', type=int, default=1000, help='Size of the \
                    region to look up for exact matching in the original \
                    assembly. Default: 1000bp')

args = parser.parse_args()

# CUSTOM FUNCTIONS


def eprint(*args, **kwargs):
    # This function prints its arguments to the standard error
    print(*args, file=sys.stderr, **kwargs)


def findall_str(sub, string):
    """
    Finds All occurrences of sub in string and returns the index of each.
    sub : Query substring
    string : String in which the lookup is performed
    returns : an iterator with all indexes of the occurrences
    """
    def next_index(length):
        index = 0 - length
        while True:
            # At every iteration, start looking after last occurrence
            index = string.find(sub, index + length)
            yield index
    # Build generator until hitting -1 (pattern not found)
    return iter(next_index(len(sub)).next, -1)

# Handling query position, whether it is from CL arg or stdin


try:
    pos1 = args.pos1.read()
except AttributeError:
    pos1 = args.pos1
# Removing potential newline
pos1 = pos1.rstrip('\n')

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
        shift_sub = [0, len(chrom.seq)]
        start_sub = bp - args.region_size / 2
        end_sub = bp + args.region_size / 2
        offset = args.region_size / 2
        if args.region_size <= len(chrom.seq):
            # Region fits in chromosome, shift if needed
            region_size = args.region_size
            shift_sub[0] -= start_sub
            shift_sub[1] -= end_sub
            if shift_sub[0] > 0:
                # If start had a negative value
                start_sub = 0  # Shift start to zero
                end_sub += shift_sub[0]  # Shift end to the right accordingly
                # Offset from start to query position
                offset = bp
            if shift_sub[1] < 0:
                # If end was larger than the length of chromosome
                end_sub = len(chrom.seq)  # Shifting end to length of chrom.
                start_sub += shift_sub[1]  # Shifting start to the left
                # Offset from start to query position
                offset = args.region_size / 2 - shift_sub[1]
        else:
            # Region does not fit in chromosome, trim it
            offset = bp
            start_sub = 0
            end_sub = len(chrom.seq)
            region_size = len(chrom.seq)
            eprint("Warning: The lookup region did not fit inside the \
chromosome and has been trimmed down to {0} bp.".format(len(chrom.seq)))
        lookup_seq = chrom.seq[start_sub:end_sub]
        break

"""
# Optional: Write query region to file.
try:
    fasta_out = open("store_region.fasta",'a')
    fasta_out.write('>'+Chr+'_'+str(bp)+'\n')
    fasta_out.write(str(lookup_seq)+'\n')
    fasta_out.close()
except NameError:
    eprint("Error: Your position was not found in the matching assembly.")
    exit()
"""
# Setting up sequence with all possible combination of rev. and comp.
# Associating boolean flag with each possibility to remember if seq. was rev.
lookups = [[lookup_seq, 0],
           [lookup_seq.complement(), 0],
           [lookup_seq.reverse_complement().complement(), 1],
           [lookup_seq.reverse_complement(), 1]]

seq_match = []
mod_flag = []

for comb_id, comb_seq in enumerate(lookups):
    # Iterate over all contigs in original assembly
    for contig in SeqIO.parse(args.ref2, "fasta"):
        # coord_match = contig.seq.find(comb_seq[0])
        # coord_match = findall_str(comb_seq[0],contig.seq)
        # iterate over all iterations of lookup sequence in current contig
        for tighit in findall_str(comb_seq[0], contig.seq):
            if tighit >= 0:
                # if the region is found in the contig, record contig name
                # and bp position corresponding to the center of the region
                if comb_seq[1]:
                    # If sequence was reversed:
                    seq_match.append((str(contig.id),
                                      str(tighit + region_size - offset)))
                    if comb_id == 2:
                        contig_mod = 'reversed '
                        mod_flag.append('rev')
                    else:
                        contig_mod = 'reverse-complemented '
                        mod_flag.append('revcomp')
                else:
                    # If sequence orientation is conserved:
                    seq_match.append((str(contig.id),
                                      str(tighit + offset)))
                    if comb_id == 1:
                        contig_mod = 'complemented '
                        mod_flag.append('comp')
                    else:
                        contig_mod = ''
                        mod_flag.append('None')
                eprint('Match found in a {0}contig !'.format(contig_mod))

if len(seq_match) > 1:
    for idx, match in enumerate(seq_match):
        print(','.join(match) + ',' +
              ','.join([str(region_size), mod_flag[idx]]))
    eprint("Warning: {0} occurences of the lookup region were found in the\
 original assembly, you should increase --region_size.".format(len(seq_match)))
elif len(seq_match) == 1:
    eprint("Success: 1 corresponding coordinate found: {0}".format(','.join(seq_match[0])))
    print(','.join(seq_match[0]) + ',' +
          ','.join([str(region_size), mod_flag[0]]))
else:
    eprint("Error: Query position not found in reference.")
