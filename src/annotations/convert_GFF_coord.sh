# This script converts the coordinates in the original GFF file to coordinates
# corresponding to the new assembly. It uses the correspondence file generated
# by corresp_contigs.sh to retrieve contig coordinates in the new assembly.
# Cyril Matthey-Doret
# 02.11.2017

gff="../../data/annotations/OGS1.0_20170110.gff"
annot="../../data/annotations/blast2go.txt"
out_annot="../../data/annotations/ordered__CMD_tracks.gff"

# iterate over lines
# for each line -> grep contig in corresp file
# sed s/// line to match new coords (rev, comp, shift)
# write edited line to new GFF
