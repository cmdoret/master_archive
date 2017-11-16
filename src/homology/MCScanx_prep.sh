#!/user/bin/env bash

# This script allows to quickly setup input files for MCScanX. It requires a set
# of gene coordinates in GFF format and a genome assembly in fasta format. If
# gene coordinates need to be converted from another assembly, the old assembly
# can be specified as well. Both assemblies must contain the same sequence, but
# contigs can be merged, reordered or reverse complemented.
# 12.11.2017
# Cyril Matthey-Doret


# Help message
function usage () {
   cat <<EOF
Usage: `basename $0` -g genes -o output_folder -r ref [-p prev_ref] [-c corresp] [-l] [-h]
   -g   gtf file with genes coordinates
   -o   output folder for MCScanX files
   -r   reference genome
   -c   contig correspondance file, if gff coordinates need to be converted
   -p   previous reference genome, if gff coordinates need to be converted and
        contig correspondance file is not available.
   -l   local run. If specified, will not use LSF bsub command
   -h   displays this help
EOF
   exit 0
}



# Parsing CL arguments
while getopts ":g:o:r:p:c:lh" opt; do
   case $opt in
   g )  GFF=${OPTARG} ;;
   o )  OUT_F=${OPTARG} ;;
   r )  REF=${OPTARG};;
   p )  PREV=${OPTARG};;
   c )  CORRESP_GFF=${OPTARG};;
   l )  local=yes;;
   h )  usage ;;
   \?)  usage ;;
   esac
done

# Testing if mandatory have been provided
if [ "x" == "x$REF" ] || [ "x" == "x$GFF" ] || \
   [ "x" == "x$OUT_F" ];
then
  echo "Error: At least a reference genome, input GFF file and output folder \
  are required."
  usage
  exit 0
fi

#1: Extract records for features of interest from gff file
echo -n "Extracting transcripts coordinates from the gff file..."
MC_GFF="$OUT_F/MCScanX_genes.gff"
grep "transcript" $GFF > $MC_GFF
echo "...transcripts extracted"
#2: If necessary, convert coordinates of GFF file to new assembly
if [ "x" == "x$PREV" ] || [ "x" == "x$CORRESP" ];
then
  echo -n "Converting transcripts coordinates to new assembly..."
  OUT_GFF="${MC_GFF%.*}_conv.gff"

  bash src/convert_coord/convert_GFF_coord.sh -i "$MC_GFF" \
                                              -o "$OUT_GFF" \
                                              ${PREV:+-O "$PREV"} \
                                              ${REF:+-N "$REF"} \
                                              ${CORRESP_GFF:+-c "$CORRESP_GFF"} \
                                              ${local:+-l}
 echo "coordinates converted !"
else
  OUT_GFF=$MC_GFF
fi

#3: convert GFF to BED format and extract sequence
echo "Converting GFF to BED."
#MC_BED="${GFF_OUT%.*}.bed"
#awk 'BEGIN{OFS="\t"}
#     {id_match=match($12,/[^;]*/)
#      id=substr($12,RSTART,RLENGTH)
#      print $1,$4,$5,id,0,$7}' $MC_GFF > $MC_BED

#module add UHTS/Analysis/BEDTools/2.26.0
#MC_SEQ="$OUT_F/MCScanX_seq.fasta"
echo "Extracting nucleotide sequences from reference."
#bedtools getfasta -fi $REF -bed $MC_BED > $MC_SEQ
#module rm UHTS/Analysis/BEDTools/2.26.0;

#4: build blast database from sequences and all vs all blast
