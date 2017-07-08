#!bin/bash

# This script parses the output of cstacks into a table to determine the best value for minimum coverage.
# Cyril Matthey-Doret
# 11.04.2017

cd "$(dirname "$0")"
cdir=../../data/cstacks/
echo 'mismatch,mean_loci,mean_alleles' > cstats.csv
for d in $cdir/*;
do
    f=$d/*allel*

    # Number of loci
    L=$(zcat $f | cut  -f3 | sort | uniq | wc -l);
    # Total number of alleles
    A=$(zcat $f | wc -l);

    echo ${d##*mm-}','$L','$A >> cstats.csv;
done

cp cstats.csv ../../reports/lab_book/
