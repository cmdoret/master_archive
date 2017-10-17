# Best not changing DAT
DAT=./data
LOG=$(DAT)/logs
PROC=$(DAT)/processed


## Mapping parameters for bwa
BWA-SRC=src/mapping/desktop_versions/bwa_nix.sh
BWA-LSF=src/mapping/bwa_lsf.sh
# aln: mismatches (MM)
MM=4

## STACKS parameters
# pstacks: minimum coverage (M)
P-SRC=src/stacks_pipeline/desktop_versions/pstacks_nix.sh
P-LSF=src/stacks_pipeline/pstacks_lsf.sh
M=3

# ustacks: number of mismatches allowed (MM) and minimum coverage (M)
U-SRC=src/stacks_pipeline/desktop_versions/ustacks_nix.sh
U-LSF=src/stacks_pipeline/ustacks_lsf.sh

# cstacks: loci mismatches (LM)
C-SRC=src/stacks_pipeline/desktop_versions/cstacks_nix.sh
C-LSF=src/stacks_pipeline/cstacks_lsf.sh
LM=3

# sstacks:
S-SRC=src/stacks_pipeline/desktop_versions/sstacks_nix.sh
S-LSF=src/stacks_pipeline/sstacks_lsf.sh

# populations: Proportion of individuals with locus (R), Minimum stack depth (D)
GR-SRC=src/stacks_pipeline/group_sstacks.sh
POP-SRC=src/stacks_pipeline/desktop_versions/populations_nix.sh
POP-LSF=src/stacks_pipeline/populations_lsf.sh
# Boolean variable: group families in 1 populations run or split them
GRFAM=T
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
BLACK=$(DAT)/ploidy/blacklist.tsv
# Note: threshold for ploidy threshold in haplo_males for my samples are:
# 0.77 if running pop. per fam.; 0.90 if running pop. on all samples

# Association mapping:
ASSOC-SRC=src/assoc_mapping/
# Number of CSD loci considered
NCSD=2


MAP=$(DAT)/mapped/aln-$(MM)
PSTACK=$(DAT)/pstacks/covmin-$(M)
USTACK=$(DAT)/ustacks/cov-$(M)_mm-$(MM)
CSTACK=$(DAT)/cstacks/mm-$(LM)
SSTACK=$(DAT)/sstacks
POP=$(DAT)/populations/grouped_d-$(D)_r-$(R)
ASSOC=$(DAT)/assoc_mapping
CENTRO=$(ASSOC)/centro
