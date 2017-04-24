MAIN=/scratch/beegfs/monthly/cmatthey
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
LM=1

# sstacks:
S-SRC=src/stacks_pipeline/multi_sstacks.sh

MAP=$(DAT)/mapped/$(ALG)-$(MM)-$(K)-$(W)
PSTACK=$(DAT)/pstacks/covmin-$(M)
CSTACK=$(DAT)/cstacks/mm-$(LM)
SSTACK=$(DAT)/sstacks
