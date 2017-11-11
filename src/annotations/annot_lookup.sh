#!/bin/bash
# This is a simple script to look at annotations nearby my GWAS hits.
# Cyril Matthey-Doret
# 20.10.2017

gff="../../data/annotations/ordered_CMD_tracks.gff"
annot="../../data/annotations/blast2go.txt"
hits='../../data/assoc_mapping/case_control/case_control_hits.tsv'
out_annot='../../data/annotations/csd_gwas_hit_annot.tsv'

#echo "Processing coordinates: chr4,1798" > casper_query.log
#python2 chr2contig.py $order_ref $orig_ref --pos1 chr4,1798  --region_size 100 2>> casper_query.log |
#python2 pos2annot.py $gff $annot --range_size 500000 > casper_query.txt 2>> casper_query.log

# reinitialize output file
rm "store_region.fasta"
echo -e "chrom\tstart\tend\tID\tGO\tterm" > $out_annot
# Convert ordered assembly coordinates to original assembly and look up
# annotations around the coordinate
while read line
do
  echo "-------------------------------"
  echo -n "Processing coordinates: "
  echo $line |
  awk 'BEGIN{RS="\t";OFS=",";}{print $2,$3}' |
  tee /dev/stderr |
  python2 ../convert_coord/pos2annot.py $gff $annot --range_size 200000 >> $out_annot
done < <(tail -n +2 $hits)
