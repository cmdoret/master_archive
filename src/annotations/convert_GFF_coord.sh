# This script converts the coordinates in the original GFF file to coordinates
# corresponding to the new assembly. It uses the correspondence file generated
# by corresp_contigs.sh to retrieve contig coordinates in the new assembly.
# Cyril Matthey-Doret
# 02.11.2017

#!/bin/bash
#BSUB -J conv_GFF
#BSUB -q normal
#BSUB -e ../../data/logs/gff_corresp.err
#BSUB -o ../../data/logs/gff_corresp.out
#BSUB -M 16000000
#BSUB -R "rusage[mem=16000]"
#BSUB -n 28
#BSUB -R "span[ptile=28]"

MAX_PROC=24
gff="../../data/annotations/OGS1.0_20170110.gff"
annot="../../data/annotations/blast2go.txt"
out_gff="../../data/annotations/ordered_CMD_tracks.gff"
corresp_gff='../../data/annotations/corresp_gff.csv'

echo -n "" > $out_gff
# iterate over lines
while read line
do
  while [ `jobs | wc -l` -ge $MAX_PROC ]
  do
    sleep 1
  done
  ( track=( $line )
  corresp=( $(grep "^${track[0]}" $corresp_gff | sed 's/,/ /g') )
  track[0]=${corresp[1]}
  # Shift start and end if necessary and flip,
  # depending if contig was reversed or not
  if [[ ${corresp[4]} == *"rev"* ]]
  then
    # Reversed: corresp[2] will match the end of the contig
    # start -> contig_end - track_end, end -> contig_end - track_start
    start=$((${corresp[2]}-${track[4]}))
    end=$((${corresp[2]}-${track[3]}))
    track[3]=$start
    track[4]=$end
  else
    # not reversed: start shifted, end -> start + contig size
    let "track[3] += ${corresp[2]}"
    let "track[4] += ${corresp[2]}"
  fi
  # If complementary and not strand agnostic -> complement track
  if [[ ${corresp[4]} == *"comp"* && ${track[6]} != "." ]]
  then
    # Reverse strand
    if [[ ${track[6]} == "+" ]]
    then
      track[6]="-"
    else
      track[6]="+"
    fi
  fi    # redirect line to output gff (line >> file)
  echo "${track[@]}" | tr ' ' \\t >> $out_gff ) &
done < $gff
wait

sort -k1,1 -k4,4n -o $out_gff $out_gff
