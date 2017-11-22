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

# If on cluster, use bsub to submit jobs, otherwise run directly
if [ -z ${local+x} ];then run_fun="bsub -I -tty";else run_fun="";fi

#1: Extract records for features of interest from gff file
echo -n "Extracting transcripts coordinates from the gff file..."
MC_GFF="$OUT_F/MCScanX_genes.gff"
awk '$3 ~ "transcript" {print $0}' "$GFF" > "$MC_GFF"
echo "transcripts extracted !"

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
# Extracting only transcripts that mapped to ordered corresp_contigs
# (i.e. in chromosomes) and removing double quotes from identifiers
awk 'BEGIN{FS="\t"} $1 ~ "chr" {print $0}' "$OUT_GFF" | \
  sed 's/"//g' > "$OUT_GFF.tmp" && mv "$OUT_GFF.tmp" "$OUT_GFF"

#3: convert GFF to BED format and extract sequence
echo "Converting GFF to BED."
OUT_BED="${OUT_GFF%.*}.bed"
awk 'BEGIN{OFS="\t"}
     {id_match=match($12,/[^;]*/)
      id=substr($12,RSTART,RLENGTH)
      print $1,$4,$5,id,0,$7}' $OUT_GFF > $OUT_BED

awk 'BEGIN{OFS="\t"}{print $1,$4,$2,$3}' "$OUT_BED" > "$OUT_GFF"
OUT_SEQ="$OUT_BED.fasta"

# If local is declared, expands to "x", otherwise, to nothing
if [ -z ${local+x} ]
then
  module add UHTS/Analysis/BEDTools/2.22.1
  module add Blast/ncbi-blast/2.6.0+
  export PATH="~/scratch/softwares/MCScanX/:$PATH"
else
  export PATH="~/Public/ncbi-blast-2.7.1+/bin:$PATH"
  export PATH="~/Public/bedtools2/bin:$PATH"
  export PATH="~/Public/MCScanX/:$PATH"
fi
echo "Extracting gene sequences from reference genome"
eval $run_fun "bedtools getfasta -fi $REF -bed $OUT_BED -fo $OUT_SEQ -name"

#4: build blast database from sequences and all vs all blast
eval $run_fun "makeblastdb -in $OUT_SEQ -dbtype nucl"

echo "Blasting transcriptome against itself."
eval $run_fun "blastn -query $OUT_SEQ \
                      -db $OUT_SEQ \
                      -outfmt 6 \
                      -max_target_seqs 5 \
                      -out $OUT_SEQ.blast"

# Renaming chromosomes according to MCScanX conventions
MC_IN="$OUT_F/MCScanX_in"
sed 's/^chr\([0-9]*\)/lf\1/' "$OUT_GFF" > "$MC_IN.gff"
sed 's/chr\([0-9]*\)/lf\1/g' "$OUT_SEQ.blast" > "$MC_IN.blast"

echo "All input files are ready:"
echo "  - BLAST output: $MC_IN.blast"
echo "  - GFF file: $MC_IN.gff"

echo "Running MCScanX"
eval $run_fun "MCScanX $MC_IN"

echo "Generating graphics control file for circle plotter"
# 800 pixels, displaying all chromosomes in the input GFF file
echo "800" > $OUT_F/circle.ctl
cut -f1 "$MC_IN.gff" | uniq | paste -s -d, - >> $OUT_F/circle.ctl
