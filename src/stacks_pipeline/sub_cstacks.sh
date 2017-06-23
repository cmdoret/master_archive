#!/bin/bash

#BSUB -o %J_STDOUT.log
#BSUB -e %J_STDERR.log
#BSUB -J Cstacks_${j}
#BSUB -n 3
#BSUB -R "span[ptile=3]"
#BSUB -R "rusage[mem=30000]"
#BSUB -M 30000000
#BSUB -q normal

## This script will build a catalogue of loci using cstacks, given a set of pstacks output.
## It will exclude samples with very low number of reads from the catalogue (less than 10% of
## the mean number of reads across samples)
# Cyril Matthey-Doret
# 22.04.2017

M=3
wd=/scratch/beegfs/monthly/cmatthey/data ## the directory containing the pstacks files
MM=1  # Mismatches allowed between samples when building loci
declare -i n=0 tot=0  # n is the number of samples, tot will store total number of radtags

for f in $wd/pstacks/covmin-$M/*tags*;  # Iterating over all samples
do
    tot+=$(zcat $f | wc -l);  # adding number of radtags to total
    n+=1;  # incrementing number of samples by one at each iteration
done;


samp=""  # initiating list of "good" samples
for i in $wd/pstacks/covmin-$M/*tags*;  # Iterating over samples (again)
do
    if [ "$(zcat $i | wc -l)" -gt $(($tot/($n*10))) ];
    then
        samp+="-s ${i%%.tags*} "  # only samples containing more than 10% of (arithmetic) mean radtags are used
    fi;
done;

module add UHTS/Analysis/stacks/1.30;

cstacks -o $wd/cstacks/mm-$MM -n $MM -p 3 $samp;
