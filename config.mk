# Best not changing DAT
DAT=./data
LOG=$(DAT)/logs
PROC=$(DAT)/processed


## Mapping parameters for bwa
BWA-SRC=src/mapping/bwa_nix.sh
BWA-LSF=src/mapping/bwa_lsf.sh
# aln: mismatches (MM)
MM=4

## STACKS parameters
# pstacks: minimum coverage (M)
P-SRC=src/stacks_pipeline/pstacks_nix.sh
P-LSF=src/stacks_pipeline/pstacks_lsf.sh
M=3
# cstacks: loci mismatches (LM)
C-SRC=src/stacks_pipeline/cstacks_nix.sh
C-LSF=src/stacks_pipeline/cstacks_lsf.sh
LM=3

# sstacks:
S-SRC=src/stacks_pipeline/sstacks_nix.sh
S-LSF=src/stacks_pipeline/sstacks_lsf.sh

# populations: Proportion of individuals with locus (R), Minimum stack depth (D)
GR-SRC=src/stacks_pipeline/group_sstacks.sh
POP-SRC=src/stacks_pipeline/populations_nix.sh
POP-LSF=src/stacks_pipeline/populations_lsf.sh
# Boolean variable: group families in 1 populations run or split them
GRFAM=T
R=80
D=5

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


MAP=$(DAT)/mapped/aln-$(MM)
PSTACK=$(DAT)/pstacks/covmin-$(M)
CSTACK=$(DAT)/cstacks/mm-$(LM)
SSTACK=$(DAT)/sstacks
POP=$(DAT)/populations/grouped_d-$(D)_r-$(R)
ASSOC=$(DAT)/assoc_mapping
CENTRO=$(ASSOC)/centro
