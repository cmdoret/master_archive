#!/bin/bash 
#BSUB -J POP
#BSUB -q normal
#BSUB -e data/logs/populations/POP_STDERR.err
#BSUB -o data/logs/populations/POP_STDOUT.out
#BSUB -M 20000000
#BSUB -R "rusage[mem=20000]"
#BSUB -n 12
#BSUB -R "span[ptile=12]"

# This script uses the populations component from the STACKS suite to compute FST across the genome of individuals. 
# It requires stacks files for each individual and a population map.
# Cyril Matthey-Doret
# 24.04.2017

od=data/populations/r-75
R=0.75

module add UHTS/Analysis/stacks/1.30;

populations -P data/sstacks/ -M data/popmap -k -b 0 -r $R -f p_value -t 12 --verbose --fstats --bootstrap --vcf

module add UHTS/Analysis/stacks/1.30;
