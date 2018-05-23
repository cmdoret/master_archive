#!bin/bash

# This script parses the output of cstacks into a table to determine the best value for minimum coverage.
# Cyril Matthey-Doret
# 11.04.2017
cd $(dirname "$0")
cdir="../../data/cstacks/"

echo 'mismatch,N_loci,N_alleles' > cstats.csv
for d in $cdir/*;
do
    f=$d/*allel*

    # Number of loci
    L=$(gunzip -c $f | cut  -f3 | sort | uniq | wc -l);
    # Total number of alleles
    A=$(gunzip -c $f | wc -l);

    echo ${d##*mm-}','$L','$A >> cstats.csv;
done

cp cstats.csv ../../docs/lab_book/cstats.csv
