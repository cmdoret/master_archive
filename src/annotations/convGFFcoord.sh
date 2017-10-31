# The purpose of this script is to convert all coordinates in a GFF file to
# matching coordinates in a reordered assembly.
# Cyril Matthey-Doret
# 31.10.2017

# For each row in the GFF, look up matching region in different assemblies
# EDIT: Would be way too slow: 192k records -> instead run once per contig and
# record contig transformation, e.g: chr3, revcomp, shift:+1059bp

# generate new GFF with matched regions

# Best possible solution: reverse transformation from jens linkage map, but no
# info on complment / reverse
