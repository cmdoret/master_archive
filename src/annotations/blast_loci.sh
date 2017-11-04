#!/bin/bash
# This script takes a list of genomic positions (columns 2 and 3 need to be
# chromosome and position in chromosome, respectively) and blast the regions
# surrounding these positions against some database, storing the blast output
# in a separate file.
# 11.07.2017
# Cyril Matthey-Doret

loci=$(cut -f2,3  -d, $1 |tail -n +2)  # Selecting chromosome and bp columns without header
genome=$2  # path of genome fasta file
len_query=1000  # desired size of sequences to blast
out_path=../../data/assoc_mapping/blast_files
mkdir -p $out_path
rm -rf $out_path/query.fasta
module add Blast/ncbi-blast/2.6.0+

for pos in $loci;  # iterating over queries
do
    chr=${pos%%,*}  # chromosome containing query
    bp=${pos##*,}  # basepair position within chromosome
    sed -n -e "/>$chr/,/>/ p" $genome | sed -e '1d;$d' | tr -d '\n' > $out_path/seq
    # extracting chromosome sequence and removing header and newlines
    start=$(($bp-$len_query/2))  # start position of query in chromosome
    query=$(awk "{print substr(\$0,$start,$len_query)}" $out_path/seq)
    echo ">"$pos >> $out_path/query.fasta
    echo $query >> $out_path/query.fasta
    # extracting query using indexes
done
bsub -K -o "blast_log.out" -e "blast_log.err" -M 20000000 -J "blast" "blastn -task dc-megablast \
-query $out_path/query.fasta -db refseq -outfmt 6 -max_target_seqs 1 -out blast.out"
module rm Blast/blast/2.2.26
