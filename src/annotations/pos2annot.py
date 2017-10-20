"""
This script takes a genomic position as input (chr BP) and a gff files based
on the same assembly and returns either the annotations within a given region
or the top N closest annotations.
Cyril Matthey-Doret
20.10.2017
"""

from __future__ import print_function  # Use python3 print function
import sys  # Will allow to print messages to stderr
import argparse  # Parses command line arguments
import pybedtools as bed  # Python wrapper for bedtools

# Parsing command line arguments

parser = argparse.ArgumentParser(description="Takes a genomic position in the \
                                 form 'Chr BP' where Chr is the name of the \
                                 contig and BP is an integer representing the \
                                 basepair position within the contig. Returns \
                                 the annotations in the neighborhood of that \
                                 position. A GFF files with annotations must \
                                 be provided.")

parser.add_argument('Chr',type=str,help="The name of the contig where the \
                    position is.")
parser.add_argument('bp',type=int,help="An integer number representing the \
                    genomic position within the contig, in basepairs.")
parser.add_argument('annot', type=str, help="The path to the GFF file \
                    containing the annotations.")
parser.add_argument('--method', type=str, choices=['range','top'],
                    default='range', help='Method used to find annotations, \
                    can be either "range", in which case all annotations \
                    within a basepair range (default: 10kb) will be picked up, \
                    or "top", in which case the top N (default: 10) nearest \
                    annotations will be used.')
parser.add_argument('--range_size', type=int, default=10000, help='In case the \
                    method used it "range", an integer defining the range in \
                    which to look for annotations around the input position.' )
parser.add_argument('--top_count', type=int, default=10, help='In case the \
                    method used is "top", the number of closest neighbour \
                    annotations that should be returned.')

args = parser.parse_args()


# Reading GFF file
