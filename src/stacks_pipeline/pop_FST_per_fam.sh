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

module add UHTS/Analysis/stacks/1.30;

for fam in $(cut -f3 data/individuals | sort | tail -n +2 | uniq))
# All families in dataset (excluding header with tail)
do
    echo "Populations running on family $fam..."
    mkdir -p $od/$fam  # Prepare one output folder per family
    populations -P data/sstacks/$fam -M data/popmap -p 2 -m $D -b 0 -r $R -f p_value -t 3 --verbose --fstats --vcf --max_obs_het 0.9
    mv data/sstacks/$fam/batch* $od/$fam/
    # Moving all populations output file from sstacks family folder to populations family folder
done    


module add UHTS/Analysis/stacks/1.30;
