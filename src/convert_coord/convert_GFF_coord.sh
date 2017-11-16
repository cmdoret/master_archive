#!/bin/bash
# This script converts the coordinates in the original GFF file to coordinates
# corresponding to the new assembly. It uses the correspondence file generated
# by corresp_contigs.sh to retrieve contig coordinates in the new assembly.
# Cyril Matthey-Doret
# 02.11.2017

# Help message
function usage () {
   cat <<EOF
Usage: `basename $0` -i input_gff -o output_gff [-c corresp_file] [-O old] [-N new] [-l] [-h]
   -i   gff file to be converted
   -o   output converted gff file
   -c   csv file for correspondance between contigs
   -O   old reference, only needed if -c is not specified
   -N   new reference, only needed if -c is not specified
   -l   local run. If specified, will not use LSF bsub command
   -h   displays this help
EOF
   exit 0
}

# Parsing CL arguments
while getopts ":i:o:c:O:N:lh" opt; do
   case $opt in

   i )  GFF=${OPTARG} ;;
   o )  OUT_GFF=${OPTARG} ;;
   c )  CORRESP_GFF=${OPTARG};;
   O )  OLD_REF=${OPTARG};;
   N )  NEW_REF=${OPTARG};;
   l )  local=yes;;
   h )  usage ;;
   \?)  usage ;;
   esac
done

if [ "x" == "x$GFF" ] || [ "x" == "x$OUT_GFF" ];
then
  echo "Error: Input and output GFF files must be provided along with either a \
  contig correspondance file or both references."
  usage
  exit 0
fi

shift $(($OPTIND - 1))
# set command to be used depending if running locally or on LSF
if [ -z ${local+x} ];then run_fun="bsub -I -tty";else run_fun="bash";fi

# Generate correspondance file if not specified
if [ -z ${CORRESP_GFF+x} ];
then
  if [ -z "$OLD_REF" ] || [ -z "$NEW_REF" ];
  then
    echo "Error: If no correspondance file is provided, both the old and new \
    references are required."
    usage
    exit 0
  fi
  CORRESP_GFF="data/annotations/corresp_gff.csv"
  bash corresp_contigs.sh -O "$OLD_REF" \
                          -N "$NEW_REF" \
                          -c "$CORRESP_GFF" ${local:+-l}
else if [ ! -f "$CORRESP_GFF" ];
  then
    echo "Correspondance file invalid. Exiting"
    exit 0
  fi
fi

# Run main code with bsub or directly, depending on -l flag
# Note many variables are escaped to force remote expansion
echo "Job submitted !"

eval $run_fun <<CONV_COORD

  #!/bin/bash

  #BSUB -J conv_GFF
  #BSUB -q normal
  #BSUB -e data/logs/gff_conv.err
  #BSUB -o data/logs/gff_conv.out
  #BSUB -M 16000000
  #BSUB -R "rusage[mem=16000]"
  #BSUB -n 28
  #BSUB -R "span[ptile=28]"

  source src/misc/jobs_manager.sh
  MAX_PROC=24
  n_rec=\$(wc -l $GFF | awk '{print \$1}')
  echo -n "" > $OUT_GFF
  # iterate over lines
  while read line
  do
    while [ \$(jobs | wc -l) -ge \$MAX_PROC ]
    do
      sleep 1
    done

    prettyload \$(wc -l $OUT_GFF | awk '{print \$1}') \$n_rec
    ( track=( \$line )
    corresp=( \$(grep "^\${track[0]}" $CORRESP_GFF | sed 's/,/ /g') )
    track[0]=\${corresp[1]}
    # Shift start and end if necessary and flip,
    # depending if contig was reversed or not
    if [[ \${corresp[4]} == *"rev"* ]]
    then
      # Reversed: corresp[2] will match the end of the contig
      # start -> contig_end-track_end, end -> contig_end-track_start
      start=\$((\${corresp[2]}-\${track[4]}))
      end=\$((\${corresp[2]}-\${track[3]}))
      track[3]=\$start
      track[4]=\$end
    else
      # not reversed: start shifted, end -> start + contig size
      let "track[3] += \${corresp[2]}"
      let "track[4] += \${corresp[2]}"
    fi
    # If complementary and not strand agnostic -> complement track
    if [[ \${corresp[4]} == *"comp"* && \${track[6]} != "." ]]
    then
      # Reverse strand
      if [[ \${track[6]} == "+" ]]
      then
        track[6]="-"
      else
        track[6]="+"
      fi
    fi    # redirect line to output gff (line >> file)
    echo "\${track[@]}" | tr ' ' \\\\t >> $OUT_GFF ) &
  done < $GFF
  wait

  sort -k1,1 -k4,4n -o $OUT_GFF $OUT_GFF
CONV_COORD
