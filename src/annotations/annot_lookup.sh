# This is a simple script to look
# Cyril Matthey-Doret
# 20.10.2017

orig_ref="../../data/ref_genome/canu2_low_corrRE.contigs.fasta"
order_ref="../../data/ref_genome/ordered_genome/merged.fasta"
gff="../../data/annotations/OGS1.0_20170110.gff"
annot="../../data/annotations/blast2go.txt"
hits='../../data/assoc_mapping/case_control/case_control_hits.tsv'

# Convert ordered assembly coordinates to original assembly and look up
# annotations around the coordinate
while read line
do
  echo $line |
  awk 'BEGIN{RS="\t";OFS=",";}{print $2,$3}' | tee /dev/stderr
  python2 chr2contig.py $order_ref $orig_ref |
  python2 pos2annot.py $gff $annot
done < <(tail -n +2 $hits)
