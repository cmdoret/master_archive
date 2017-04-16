#!bin/bash

# This script parses the output of pstacks into a table to determine the best value for minimum coverage.
# Cyril Matthey-Doret
# 11.04.2017

module add R/latest

pdir=/scratch/beegfs/monthly/cmatthey/data/pstacks
echo 'min_cov,mean_loci,mean_alleles' > pstats.csv

for d in $pdir/*;
do
    declare -i L=0 A=0 C=0
    for f in $d/*.alleles.tsv.gz;
    do
        # Number of loci
        L+=$(zcat $f | cut  -f3 | uniq | wc -l);
        # Total number of alleles
        A+=$(zcat $f | wc -l);
        # Number of samples
        C+=1
    done;
    FL=$(echo $L/$C | R --vanilla --quiet | sed -n '2s/.* //p')
    FA=$(echo $A/$C | R --vanilla --quiet | sed -n '2s/.* //p')
    echo ${d##*covmin-}','$FL','$FA >> pstats.csv
done;

module rm R/latest
cp pstats.csv ../../reports/lab_book/
