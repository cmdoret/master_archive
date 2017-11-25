
# Extracting all BLAST hits between genes (non-self) and parse the results into
# A circos-compatible format.
# Cyril Matthey-Doret
# 24.11.2017

# Tabular output (fmt6) from BLAST
blast='data/homology/MCScanX/input/MCScanX_in.blast'
# BED file with all genes coordinates
bed="data/homology/MCScanX/input/MCScanX_genes_conv.bed"
# Significance threshold for E-values
Esig="1*10^-5"
# Output folder
circos_links='data/circos/blast.lf.txt'
source src/misc/jobs_manager.sh

echo -n "" > "$circos_links"
declare -i i=0
N=$(wc -l $blast | awk '{print $1}')
while read -a line
do
  # Excluding self hits
  if [ "${line[1]}" != "${line[2]}" ]
  then
    # Formatting scientific notation
    Eval=$(echo "${line[11]}" | sed 's/[eE]+\{0,1\}/*10^/g')
    # Testing if E-value is below significance threshold
    if [ $(echo "$Eval  < $Esig" | bc -l) ]
    then
      grep "${line[1]}" $bed |
        awk '{printf "%s\t%s\t%s\t",$1,$2,$3}' >> $circos_links
      grep "${line[2]}" $bed |
        awk '{print $1,$2,$3}' >> $circos_links

    fi
  fi
  prettyload $i $N
  ((i++))
done < $blast
