# Builds a list of contig containing GWAS hits with p-value below threshold.
# Arguments:
#  1: list of GWAS hits
#  2: significance threshold
#  3: correspondance list of contigs coordinates vs chromosome
#  4: output folder
# Cyril Matthey-Doret
# 14.02.2018


GWAS_HIT=$1
SIGPOW=$2
TIG=$3
OUT=$4

# Coordinates of contig boundaries in chromosomes
tail -n +2 $TIG |
  awk 'BEGIN{
         FS=",";OFS=" "}
       $2 ~ "chr" {
         if ( $5 ~ "rev" ){print $2,$3-$4,$3}
       else {print $2,$3,$3+$4}
       }' |
  sort -k1,1 -k2,2n > "$OUT/tig_boundaries.tsv"


# Iterate over GWAS hits
while read line
do
  # is the p-value significant ?
  pval=$(echo $line | awk '{print $NF}' | sed 's/[eE]+\{0,1\}/*10^/g')
  if (( $(echo "$pval  > $SIGPOW" | bc -l) ))
  then
    # print contig containing significant hits
    awk -v chrom="$(cut -d$' ' -f2 <(echo $line))" \
        -v pos="$(cut -d$' ' -f3 <(echo $line))" \
        '$1 == chrom {
          if ( $2 <= pos && $3 >= pos ) {print $0}}' "$OUT/tig_boundaries.tsv"
  fi
done < <( tail -n +2 "$GWAS_HIT") |
  uniq > "$OUT/CSD_tig.tsv"
