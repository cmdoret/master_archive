#!/bin/bash 
#BSUB -J POP
#BSUB -q normal
#BSUB -e data/logs/populations/POP_STDERR.err
#BSUB -o data/logs/populations/POP_STDOUT.out
#BSUB -M 20000000
#BSUB -n 3

# This script uses the populations component from the STACKS suite to compute FST across the genome of individuals. 
# It requires stacks files for each individual and a population map.
# Cyril Matthey-Doret
# 24.04.2017

od=data/populations/r-75
R=0.75

module add UHTS/Analysis/stacks/1.30;

populations -P data/sstacks/ -M data/popmap -b 1 -k -r $R -f p_value -o $od -t 3

module add UHTS/Analysis/stacks/1.30;