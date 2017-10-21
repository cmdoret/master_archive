#!/bin/bash
# This is a simple script to look at annotations nearby my GWAS hits.
# Cyril Matthey-Doret
# 20.10.2017

orig_ref="../../data/ref_genome/canu2_low_corrRE.contigs.fasta"
order_ref="../../data/ref_genome/ordered_genome/merged.fasta"
gff="../../data/annotations/OGS1.0_20170110.gff"
annot="../../data/annotations/blast2go.txt"
hits='../../data/assoc_mapping/case_control/case_control_hits.tsv'

#echo "Processing coordinates: chr4,1798" > casper_query.log
#python2 chr2contig.py $order_ref $orig_ref --pos1 chr4,1798  --region_size 100 2>> casper_query.log |
#python2 pos2annot.py $gff $annot --range_size 500000 > casper_query.txt 2>> casper_query.log

# Convert ordered assembly coordinates to original assembly and look up
# annotations around the coordinate
while read line
do
  echo "-------------------------------"
  echo -n "Processing coordinates: "
  echo $line |
  awk 'BEGIN{RS="\t";OFS=",";}{print $2,$3}' |
  tee /dev/stderr |
  python2 chr2contig.py $order_ref $orig_ref |
  python2 pos2annot.py $gff $annot > $(echo $line | awk 'BEGIN{RS="\t";OFS="_";}{print $2,$3}')
done < <(tail -n +2 $hits)
