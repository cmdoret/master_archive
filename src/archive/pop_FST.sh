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
# It requires stacks files for each individual and a population map.
# Cyril Matthey-Doret
# 24.04.2017

sstacks=data/sstacks
od=./data/populations/d-5_r-75
R=0.75
D=5

module add UHTS/Analysis/stacks/1.46;

populations -P data/sstacks/ -M data/popmap -p 2 -m $D -b 0 -r $R -f p_value -t 3 --verbose --fstats --vcf --max_obs_het 0.9
mv $sstacks/batch* $od
