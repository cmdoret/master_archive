#!/bin/bash
# extracting CSD contig, masking coding regions and blasting to find noncoding
# homologies
# Cyril Matthey-Doret
# 25.11.2017

export PATH="~/Public/ncbi-blast-2.7.1+/bin/:$PATH"
export PATH="~/Public/bedtools2/bin/:$PATH"


bedtools maskfasta -fi 'data/ref_genome/ordered_genome/merged.fasta' \
                   -bed 'data/homology/MCScanX/input/MCScanX_genes_conv.bed' \
                   -fo 'data/homology/ref_masked_coding.fasta'

awk 'BEGIN {RS=">"} /chr/ {print ">"$0}' 'data/homology/ref_masked_coding.fasta' |
  sed 's/chr/lf/g' > 'data/homology/sub_ref_masked.fasta'

#sed $'s/ /\t/g' "data/circos/csd_tig.lf.txt" |
#  sed 's/lf/chr/' > 'data/homology/csd_tig.bed'

#bedtools getfasta -fi 'data/homology/ref_masked_coding.fasta' \
#                  -bed 'data/homology/csd_tig.bed' \
#                  -fo 'data/homology/csd_noncoding.fasta'


makeblastdb -in 'data/homology/sub_ref_masked.fasta' -dbtype nucl
blastn -query 'data/homology/sub_ref_masked.fasta' \
       -db 'data/homology/sub_ref_masked.fasta' \
       -outfmt 6 -max_target_seqs 5 \
       -out 'data/homology/noncoding.blast'

awk '$1 != $2 {print $0}' "data/homology/noncoding.blast" > "data/homology/noncoding_inter.blast"
