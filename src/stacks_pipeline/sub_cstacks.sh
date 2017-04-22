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

wd=/scratch/beegfs/monthly/cmatthey/data ## the directory containing the pstacks files
MM=1
declare -i n=0 tot=0

for f in $wd/pstacks/covmin-4/*tags*;
do
    tot+=$(zcat $f | wc -l);
    n+=1;
done;


samp=""
for i in $wd/pstacks/covmin-4/*tags*;
do
    if [ "$(zcat $i | wc -l)" -gt $(($tot/($n*10))) ];
    then
        samp+="-s ${i%%.tags*} "
    fi;
done;

module add UHTS/Analysis/stacks/1.30;

cstacks -o $wd/cstacks/mm-$MM -n $MM -p 3 $samp;
