# This script allows to get a correspondance file of the coordinates of all
# contigs between 2 references. Both references must contain the same sequence,
# but contigs can be reorded, reverse-complemented and merged in the new
# reference.
# Cyril Matthey-Doret
# 31.10.2017

# Help message
function usage () {
   cat <<EOF
Usage: `basename $0` -O old_ref -N new_ref -c corresp [-l] [-h]
   -O   old reference
   -N   new reference
   -c   output correspondance file
   -l   local run. If specified, will not use LSF bsub command
   -h   displays this help
EOF
   exit 0
}

# Parsing CL arguments
while getopts ":O:N:lc:h" opt; do
   case $opt in
   O )  OLD_REF=${OPTARG} ;;
   N )  NEW_REF=${OPTARG} ;;
   c )  CORRESP_GFF=${OPTARG};;
   l )  local=yes;;
   h )  usage ;;
   \?)  usage ;;
   esac
done


if [ -z ${local+x} ];then run_fun="bsub";else run_fun="bash";fi

eval "$run_fun" <<CORR
#!/bin/bash
#BSUB -J corresp_GFF
#BSUB -q normal
#BSUB -e data/logs/gff_corresp.err
#BSUB -o data/logs/gff_corresp.out
#BSUB -M 32000000
#BSUB -R "rusage[mem=32000]"
#BSUB -n 36
#BSUB -R "span[ptile=36]"

MAX_PROC=32
orig_ref="$OLD_REF"
order_ref="$NEW_REF"
corresp_gff="$CORRESP_GFF"


# Get start position of all contigs in the ordered assembly and whether
# They have been reversed and/or complemented. Store infos in a file.

echo 'tig,chr,start,region_size,transform' > $corresp_gff


grep ">" $orig_ref | while read tig
do

  while [ `jobs | wc -l` -ge $MAX_PROC ]
  do
    sleep 1
  done
    # Setting region size longer than whole genome -> resized to contig
    ( tname=$(echo $tig | sed -n 's/>\([^ ]*\) .*/\1,/p' | tr -d '\n')
    new_coord=$(python2 chr2contig.py --pos1 "${tname}0" \
                          --region_size 1000000000 \
                          $orig_ref \
                          $order_ref)
    echo ${tname}${new_coord} >> $corresp_gff ) &
done
wait
CORR
