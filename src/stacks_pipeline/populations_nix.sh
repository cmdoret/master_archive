#!/bin/bash

# This script uses the populations component from the STACKS suite to compute various populations statistics the genome of individuals.
# It needs a popmap file as specified in the STACKS documentation: catchenlab.life.illinois.edu/stacks/comp/populations.php
# If group is set to T it requires stacks files for each individual in the same folder.
# Otherwise, stacks files need to be in separate folders for each family.
# Cyril Matthey-Doret
# 24.04.2017

while [[ "$#" > 1 ]]; do case $1 in
    # Path to Sstacks files
    --sst) sst="$2";;
    # Path to "threshold" file with ploidy information
    --thresh) cst="$2";;
    # Grouped or per-family populations run
    --group) group="$2";;
    # Minimum locus depth required
    --md) D="$2";;
    # Minimum proportion of individuals with locus
    --r) R="$2";;
    # Output populations folder
    --out) out_dir="$2";;
    # Location of log files
    --log) logs="$2";;
    *) break;;
  esac; shift; shift
done

if [ "$group" = "T" ]
then  # Families are grouped into a single populations run

    if [ -f $thresh ]  # If ploidy information already available
    then
        haplo=$(awk '$11 ~ "H" {print $1}' $thresh)  # list haploid males
        for indv in $haplo;
        do
            rm $sst/$indv*
            # Remove all files related to haploids to exclude from populations run
        done
        echo "Populations running on all individuals, excluding $(echo $haplo | wc -w) haploid males."
    fi
    mkdir -p $out_dir/  # Prepare one output folder per family
    populations -P $sst -M data/popmap.tsv -p 2 -m $D -b 1 -r $R -k -f p_value \
    -t 3 --verbose --fstats --vcf --max_obs_het 0.9 --genomic --renz ecoRI

    mv $sst/batch* $out_dir/
    # Moving all populations output file from sstacks family folder to populations family folder

else

    for fam in $(cut -f3 data/individuals.tsv | tail -n +2 | sort | uniq)
    # All families in dataset (excluding header with tail)
    do
        if [ -f $thresh ]  # If ploidy information already available
        then
            haplo=$(awk -v var=$fam '$11 ~ "H" && $3 ~ var {print $1}' $thresh)  # list haploid males for each family
            for indv in $haplo;
            do
                rm $sst/$fam/$indv*
                # Remove all files related to haploids to exclude from populations run
            done
        fi
        echo "Populations running on family $fam, excluding $(echo $haplo | wc -w) haploid males."
        mkdir -p $out_dir/$fam  # Prepare one output folder per family
        populations -P $sst/$fam -M data/popmap.tsv -p 2 -m $D -b 1 -r $R -k \
        -f p_value -t 3 --verbose --fstats --vcf --max_obs_het 0.9 --genomic --renz ecoRI

        mv $sst/$fam/batch* $out_dir/$fam/
        # Moving all populations output file from sstacks family folder to populations family folder

    done

fi
