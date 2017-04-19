#!/bin/bash

#BSUB -o %J_STDOUT.log
#BSUB -e %J_STDERR.log
#BSUB -J Cstacks_${j}
#BSUB -n 3
#BSUB -R "span[ptile=3]"
#BSUB -R "rusage[mem=30000]"
#BSUB -M 30000000
#BSUB -q normal

## Simple script for submitting pstacks jobs

wd=/scratch/beegfs/monthly/cmatthey/data ## the directory containing the pstacks files
MM=1



samp=""
for i in $wd/pstacks/covmin-4/*allel*;
do
    if [ "$(zcat $i | wc -l)" -gt "1" ];
    then
        samp+="-s ${i%%.alleles*} "
    fi;
done;

module add UHTS/Analysis/stacks/1.30;

cstacks -o $wd/cstacks/mm-$MM -n $MM -p 3 $samp;
