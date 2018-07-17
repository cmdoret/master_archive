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
import re
import gc
import gzip
import os
from itertools import izip_longest, islice
from time import time
# If available, use c version of pickle (much faster)
try:
    import cPickle as pickle
except ImportError:
    import pickle

###################
#  PARSE CL ARGS  #
###################


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
                    by default. Multiple positions can be given and need to \
                    be separated by either a newline or a pipe (|).')
parser.add_argument('ref1', type=str, help='Path to the fasta file of the \
                    assembly, where the position is known.')
parser.add_argument('ref2', type=str, help='Path to the fasta file of the \
                    assembly, where the position is unknown.')
parser.add_argument('--region_size', type=int, default=1000, help='Size of the \
                    region to look up for exact matching in the original \
                    assembly. Default: 1000bp')

args = parser.parse_args()

####################
# CUSTOM FUNCTIONS #
####################


class ChromSuffixArray:
    # Used to build a suffix array for speeding up
    # lookups in a chromosome
    def __init__(self, seq, build=True):
        self.seq = seq
        if build:
            self.array = self.build(self.seq)
            self.array = self.inverse_array(self.array)
        self.len = len(self.seq)

    def to_int_keys(self, l):
        """
        Convert a string to a list of integer keys
        representing the lexicographic order of keys.
        l: iterable of keys
        returns: a list with integer keys
        """
        seen = set()
        ls = []
        for e in l:
            if e not in seen:
                ls.append(e)
                seen.add(e)
        ls.sort()
        index = {v: i for i, v in enumerate(ls)}
        return [index[v] for v in l]

    def disk_load(self, path):
        gc.disable()
        disk_array = gzip.open(path, 'rb')
        self.array = pickle.load(disk_array)
        gc.enable()

    def build(self, s):
        """
        Builds the suffix array.
        s: string for which to build the array.
        O(n * log(n)^2)
        """
        n = len(s)
        k = 1
        line = self.to_int_keys(s)
        while max(line) < n - 1:
            line = self.to_int_keys(
                [a * (n + 1) + b + 1
                 for (a, b) in
                 izip_longest(line, islice(line, k, None),
                              fillvalue=-1)])
            k <<= 1
        return line

    def inverse_array(self, l):
        """
        l: suffix array to invert.
        """
        n = len(l)
        ans = [0] * n
        for i in range(n):
            ans[l[i]] = i
        return ans

    def search(self, P):
        """
        Search pattern in a suffix array using binary search.
        Search for entire string, not subsets. Can identify multiple matches.
        P: Pattern to look up (string).
        S: Suffix array to search.
        returns: list of integers indicating where P was found in the sequence,
        -1 otherwise.
        """

        query_len = len(P)
        left = 0
        right = self.len
        match_seq = []

        # Binary search
        while left < right:
            mid = (left+right) // 2
            pos = self.array[mid]
            suffix = self.seq[pos:(pos + query_len)]
            if P > suffix:
                left = mid + 1
            elif P < suffix:
                right = mid
            else:
                left = right

        # If pattern was found, add a match
        if self.seq[pos:(pos + query_len)] == P:
            match_seq.append(pos)

            # Look for identical sequence on the right
            adj_mid = mid + 1
            pos = self.array[adj_mid]
            while (self.seq[pos:(pos+query_len)] == P) & \
                  (adj_mid < (self.len-1)):
                match_seq.append(pos)
                adj_mid += 1
                pos = self.array[adj_mid]

            # Look for identical sequence on the left
            adj_mid = mid - 1
            pos = self.array[adj_mid]
            while (self.seq[pos:(pos+query_len)] == P) & (adj_mid > 0):
                match_seq.append(pos)
                adj_mid -= 1
                pos = self.array[adj_mid]
        # If pattern is not found, return -1
        else:
            match_seq.append(-1)

        return match_seq


def eprint(*args, **kwargs):
    # prints its arguments to the standard error
    print(*args, file=sys.stderr, **kwargs)


#################
# PARSE CL ARGS #
#################

# Handling query position, whether it is from CL arg or stdin
try:
    pos1 = args.pos1.read()
except AttributeError:
    pos1 = args.pos1
# Splitting positions by pipe or newline
pos1 = re.split('\n|\|', pos1)
# Removing empty elements (can arise in case of empty lines)
pos1 = [x for x in pos1 if x]

# Splitting query into chromosome and basepair
for i in range(len(pos1)):
    if len(pos1[i].split(',')):
        npos = pos1[i].split(',')
        pos1[i] = npos

        try:
            int(npos[1])
        except ValueError:
            eprint("Error: bp must be an integer number representing the \
    genomic position in basepairs within the chromosome.")
    else:
        eprint("Error: Query position must be in the form 'Chr,bp'.")

######################################
# BUILD SUFFIX ARRAY FOR EACH CONTIG #
######################################

ref1 = list(SeqIO.parse(args.ref1, "fasta"))
ref2 = list(SeqIO.parse(args.ref2, "fasta"))

tig_suf_arr = {}
refdir = os.path.dirname(args.ref2)
sufdir = os.path.join(refdir, 'suffix')

# Create directory for storing suffix arrays if not present
try:
    os.makedirs(sufdir)
# Do not raise error if the directory already exists
except OSError:
    if not os.sufdir.isdir(sufdir):
        raise

# If pickled suffix files already exist, use these  instead of building
try:
    eprint("Loading suffix array from disk: ", end='')
    # suffix array is compressed to reduce size on disk.
    # garbage collector (gc) disabled as it slows down loading.
    for tig in ref2:
        if re.search('chr.*', tig.id):
            with gzip.open(os.path.join(sufdir, tig.id), 'rb') as suffix:
                eprint(tig.id)
                tig_suf_arr[tig.id] = ChromSuffixArray(tig.seq, build=False)
                tig_suf_arr[tig.id].disk_load(suffix)
    eprint("Done !")
# If the file is absent, build it
except IOError:
    eprint("File does not exist !")
    eprint("Building suffix array...")

    # Store a suffix array for each contig in a dictionary
    for tig in ref2:
        eprint(str(tig.id) + "...")
        # Note unplaced contigs are excluded (i.e. only chromosomes)
        if re.search('chr.*', tig.id):
            tig_suf_arr[tig.id] = ChromSuffixArray(tig.seq)
            # Dump suffix array to disk for next run
            with gzip.open(os.path.join(sufdir, tig.id), 'wb') as suffix:
                pickle.dump(file=suffix, obj=tig_suf_arr[tig.id].array)
    eprint("Done !")

###########################
# PROCESS INPUT POSITIONS #
###########################

# Loop over input positions
for pos in pos1:
    Chr = pos[0]
    bp = int(pos[1])
    for chrom in ref1:
        # Iterating over chromosomes in the ordered assembly
        if chrom.id == Chr:
            # if correct chromosome, subset a region of defined size centered
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
                    end_sub += shift_sub[0]  # Shift end to the right
                    # Offset from start to query position
                    offset = bp
                if shift_sub[1] < 0:
                    # If end was larger than the length of chromosome
                    end_sub = len(chrom.seq)  # Shifting end to len. of chrom.
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

    # Setting up sequence with all possible combination of rev. and comp.
    # Associating boolean flag with each possibility if seq. was rev.
    lookups = [[lookup_seq, 0],
               [lookup_seq.complement(), 0],
               [lookup_seq.reverse_complement().complement(), 1],
               [lookup_seq.reverse_complement(), 1]]

    seq_match = []
    mod_flag = []

    for comb_id, comb_seq in enumerate(lookups):
        # Iterate over all contigs in original assembly
        for contig in ref2:
            # Do not look for match in unplaced contigs
            if not re.search('chr.*', contig.id):
                continue
            # coord_match = contig.seq.find(comb_seq[0])
            # coord_match = findall_str(comb_seq[0],contig.seq)
            # iterate over all occurences of lookup sequence in current contig
            start = time()
            for tighit in tig_suf_arr[contig.id].search(comb_seq[0]):
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
      assembly, you should increase --region_size.".format(len(seq_match)))
    elif len(seq_match) == 1:
        eprint("Success: 1 corresponding coordinate \
                found: {0}".format(','.join(seq_match[0])))
        print(','.join(seq_match[0]) + ',' +
              ','.join([str(region_size), mod_flag[0]]))
    else:
        eprint("Error: Query position not found in reference.")
