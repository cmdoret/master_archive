# Project folder
MAIN=/scratch/beegfs/monthly/cmatthey
# Best not changing DAT
DAT=./data
PROC=$(DAT)/processed

## Mapping parameters for bwa
BWA-SRC=src/mapping/bwa_script.sh
ALG=aln
# aln: mismatches (MM)
MM=4

# mem: min seed length(K) and band width (W)
K=19
W=100

## STACKS parameters
# pstacks: minimum coverage (M)
P-SRC=src/stacks_pipeline/multi_pstacks.sh
M=3
# cstacks: loci mismatches (LM)
C-SRC=src/stacks_pipeline/sub_cstacks.sh
LM=3

# sstacks:
S-SRC=src/stacks_pipeline/multi_sstacks.sh

# populations: Proportion of individuals with locus (R), Minimum stack depth (D)
GR-SRC=src/stacks_pipeline/group_sstacks.sh
GRFAM=T  # Boolean variable: group families in 1 populations run ?
# POP-SRC=src/stacks_pipeline/pop_FST.sh
POP-SRC=src/stacks_pipeline/populations.sh
R=80
D=20

# Misc: report generation
LAB=reports/lab_book
MISC=src/misc
REF=$(DAT)/ref_genome/ordered_genome/merged.fasta
REF-ANN=$(DAT)/ref_genome/ordered_genome/merged.fasta.ann

# Ploidy: exclude haplomales
VCFSUM=$(DAT)/ploidy/vcftools/summary_full.txt
THRESH=$(DAT)/ploidy/thresholds/fixed.tsv

# Association mapping:
ASSOC-SRC=src/assoc_mapping/
# Number of CSD loci considered
NCSD=2


MAP=$(DAT)/mapped/$(ALG)-$(MM)-$(K)-$(W)
PSTACK=$(DAT)/pstacks/covmin-$(M)
CSTACK=$(DAT)/cstacks/mm-$(LM)
SSTACK=$(DAT)/sstacks
POP=$(DAT)/populations/grouped_d-$(D)_r-$(R)
ASSOC=$(DAT)/assoc_mapping
CENTRO=$(ASSOC)/centro
