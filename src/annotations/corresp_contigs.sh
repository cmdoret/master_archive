# The purpose of this script is to quickly get the coordinates corresponding to
# all contigs start positions in the ordered assembly.
# Cyril Matthey-Doret
# 31.10.2017
#!/bin/bash
#BSUB -J corresp_GFF
#BSUB -q normal
#BSUB -e ../../data/logs/gff_corresp.err
#BSUB -o ../../data/logs/gff_corresp.out
#BSUB -M 32000000
#BSUB -R "rusage[mem=32000]"
#BSUB -n 36
#BSUB -R "span[ptile=36]"

MAX_PROC=32
orig_ref="../../data/ref_genome/canu2_low_corrRE.contigs.fasta"
order_ref="../../data/ref_genome/ordered_genome/merged.fasta"
corresp_gff='../../data/annotations/corresp_gff.csv'


# Get start position of all contigs in the ordered assembly and whether
# They have been reversed and/or complemented. Store infos in a file.

echo 'tig,chr,region_size,transform' > $corresp_gff


grep ">" $orig_ref | while read tig
do

  while [ `jobs | wc -l` -ge $MAX_PROC ]
  do
    sleep 1
  done

    ( tname=$(echo $tig | sed -n 's/>\([^ ]*\) .*/\1,/p' | tr -d '\n')
    new_coord=$(python2 chr2contig.py --pos1 "${tname}0" \
                          --region_size 1000000000 \
                          $orig_ref \
                          $order_ref)
    echo ${tname}${new_coord} >> $corresp_gff ) &
done
wait
