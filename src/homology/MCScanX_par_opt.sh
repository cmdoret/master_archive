#!/bin/bash
# Choosing optimal values for advanced parameters of MCScanX based on the size
# of regions of interest and number of overlapping transcripts.
# Cyril Matthey-Doret
# 17.11.2017

# Help message
function usage () {
   cat <<EOF
Usage: `basename $0` -a assoc_hits -r ref -p prev_ref -c corresp [-l] [-h]
   -a   file containing significant association mapping hits
   -r   reference genome in fasta format
   -p   previous reference in fasta format
   -c   contigs correspondance file
   -h   displays this help
EOF
   exit 0
}

# Parsing CL arguments
while getopts ":a:r:p:c:h" opt; do
   case $opt in
   a )  ASSOC=${OPTARG} ;;
   r )  REF=${OPTARG} ;;
   p )  PREV=${OPTARG};;
   c )  CORRESP=${OPTARG};;
   h )  usage ;;
   \?)  usage ;;
   esac
done

echo "Estimating optimal value for gap size parameter:"
tigs=()
while read line
do
  tig=$(echo -n "$line" | awk 'BEGIN{OFS=","}{print $2,$3}' |
    python2 src/convert_coord/chr2contig.py $REF \
       "$PREV" 2> /dev/null |
    awk 'BEGIN{FS=","}{print $1}')
    if [[ ! " ${tigs[@]} " =~ " ${tig} " ]]; then
      tigs+=( "$tig" )
      echo -ne "  1. Retrieving original contigs of GWAS hits: " "${tigs[@]}" "\r"
    fi
done < <(tail -n +2 $ASSOC)

echo -e "\n"
echo "show tigs: " "${tigs[@]}"
for tig in "${tigs[@]}"
do
  # a) Retrieve contig size in corresp_gff
  len=$(awk -v tig=$tig 'BEGIN{FS=","} $1 == tig {print $4}' $CORRESP)
  # c) count transcripts in contig
  n_genes=$(awk -v tig=$tig 'BEGIN {FS="\t";N=0}
                             $1 ~ tig && $3 == "transcript" {N++}
                             END {print N}' data/rna_seq/assembled/transcripts.gtf)
  # d) take max (or min ?) value as number of overlapping genes parameter
  echo "  2. contig stats: $tig is $len bp long and contains $n_genes transcripts"
done
