
DAT_FILES=./data
PROC=$(DAT_FILES)/processed
MAP=$(DAT_FILES)/mapped
STACK=$(DAT)/pstacks

## Mapping parameters for bwa
ALG=aln
# aln: mismatches (MM)
MM=4
# mem: min seed length (K) and band width (W)
K=19
W=100

## STACKS parameters
# pstacks: minimum coverage (M)
M=3

