# This script blasts the chromosomes containing GWAS hits against one-another.
# Cyril Matthey-Doret
# 04.11.2017

gwas='../../data/assoc_mapping/case_control/case_control_hits.tsv'
blast_files="../../data/homology/blast_files/"
mkdir -p $blast_files

module add Blast/ncbi-blast/2.6.0+;

# Extracting chromosomes containing sighificant GWAS hits
chr_hits=$(cut -f2 gwas | tail -n +2 | sort | uniq)

# Generating one fasta file for each chromosome and building local blast db
for chrom in $chr_hits
do
  awk "BEGIN{RS='>'}/$chrom/{print}" > "$blast_files/$chrom.fasta"
  makeblastdb -in "$blast_files/$chrom.fasta" -dbtype nucl
done

# Blasting chromosomes against each other (one-way)
blasted=''
for chromQ in $chr_hits
do
  for chromS in $chr_hits
  do
    # Do not blast chrom against itself of use twice the same pair
    if [[ "*${chromQ}${chromS}*" != $blasted && \
          "*${chromS}${chromQ}*" != $blasted || \
          $chromQ != $chromS ]]
     then
       blasted+="${chromQ}${chromS} "
       blastn -query "$blast_files/$chromQ.fasta" \
              -db "$blast_files/$chromS" \
              -outfmt 6 \
              -max_target_seqs 1 \
              -out "$blast_files/blast_${chromQ}_${chromS}.out"

module rm Blast/ncbi-blast/2.6.0+;

for b_out in $blast_files/*.out;
do
  
done
