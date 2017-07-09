#!/bin/bash 
#BSUB -J POP
#BSUB -q normal
#BSUB -e data/logs/populations/POP_STDERR.err
#BSUB -o data/logs/populations/POP_STDOUT.out
#BSUB -M 16000000
#BSUB -R "rusage[mem=16000]"
#BSUB -n 3
#BSUB -R "span[ptile=3]"

# This script uses the populations component from the STACKS suite to compute FST across the genome of individuals. 
# It runs the program on each family separately.
# It requires stacks files for each individual in a separate folder for each family and a population map.
# Cyril Matthey-Doret
# 24.04.2017

od=./data/populations/d-5_r-75
R=0.75
D=5
thresh=data/ploidy/thresholds/m2

module add UHTS/Analysis/stacks/1.30;

for fam in $(cut -f3 data/individuals | tail -n +2 | sort | uniq)
# All families in dataset (excluding header with tail)
do
    if [ -f $thresh ]  # If ploidy information already available
    then
        haplo=$(awk -v var=$fam '$11 ~ "H" && $3 ~ var {print $1}' $thresh)  # list haploid males for each family
        for indv in $haplo;
        do
            rm data/sstacks/$fam/$indv*
            # Remove all files related to haploids to exclude from populations run
        done
    fi
    
    echo "Populations running on family $fam, excluding $(echo $haplo | wc -w) haploid males."
    mkdir -p $od/$fam  # Prepare one output folder per family
    populations -P data/sstacks/$fam -M data/popmap -p 2 -m $D -b 0 -r $R -f p_value -t 3 --verbose --fstats --vcf --max_obs_het 0.9
    mv data/sstacks/$fam/batch* $od/$fam/
    # Moving all populations output file from sstacks family folder to populations family folder
done    

#rm $od/batch*

