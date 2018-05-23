# Best not changing DAT
DAT=./data
LOG=$(DAT)/logs
PROC=$(DAT)/processed


## Mapping parameters for bwa
BWA-RAD=src/mapping/bwa_RADseq.sh
# aln: mismatches (MM)
MM=4

## STACKS parameters
# pstacks: minimum coverage (M)
P-SRC=src/stacks_pipeline/pstacks_script.sh
M=3

# ustacks: number of mismatches allowed (MM) and minimum coverage (M)
U-SRC=src/stacks_pipeline/desktop_versions/ustacks.sh
U-LSF=src/stacks_pipeline/ustacks_lsf.sh

# cstacks: loci mismatches (LM)
C-SRC=src/stacks_pipeline/cstacks_script.sh
LM=3

# sstacks:
S-SRC=src/stacks_pipeline/sstacks_script.sh

# populations: Proportion of individuals with locus (R), Minimum stack depth (D)
GR-SRC=src/stacks_pipeline/group_sstacks.sh
POP-SRC=src/stacks_pipeline/populations_script.sh
# Boolean variable: group families in 1 populations run or split them
GRFAM=T
R=80
D=5

# Misc: report generation
LAB=docs/lab_book
MISC=src/misc
REF=$(DAT)/ref_genome/ordered_genome/merged.fasta
REF-ANN=$(DAT)/ref_genome/ordered_genome/merged.fasta.ann
SIZES=$(DAT)/ref_genome/tig.sizes.txt

# Ploidy: exclude haplomales
VCFSUM=$(DAT)/ploidy/vcftools/summary_full.txt
THRESH=$(DAT)/ploidy/thresholds/fixed.tsv
BLACK=$(DAT)/ploidy/blacklist.tsv
# Note: best value for ploidy threshold in haplo_males for our samples is 0.90
# Change it according to distribution of homozygousity in the dataset
HOM_PLOID=0.90

# Association mapping:
ASSOC-SRC=src/assoc_mapping/
# Number of CSD loci considered
NCSD=2

#RNA-seq data processing
RNA=$(DAT)/rna_seq/
BAM=$(RNA)/STAR_larvae_aligned.sort.bam.bam
RNA-SRC=src/rna_seq/process_rna.sh
OLD-REF=$(DAT)/ref_genome/canu2_low_corrRE.contigs.fasta

# MCScanX
MCSX-SRC=src/homology/MCScanX_prep.sh
MCSX-IN=$(DAT)/homology/MCScanX/input
CORRESP=$(DAT)/annotations/corresp_gff.csv

# Linkage map
LINKMAP=$(DAT)/linkage_map/


MAP=$(DAT)/mapped/aln-$(MM)
PSTACK=$(DAT)/pstacks/covmin-$(M)
USTACK=$(DAT)/ustacks/cov-$(M)_mm-$(MM)
CSTACK=$(DAT)/cstacks/mm-$(LM)
SSTACK=$(DAT)/sstacks
POP=$(DAT)/populations/grouped_d-$(D)_r-$(R)
ASSOC=$(DAT)/assoc_mapping
HITS=data/assoc_mapping/case_control/case_control_hits.tsv
CENTRO=$(ASSOC)/centro
WGS=$(DAT)/wgs_wild/
