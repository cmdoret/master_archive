!/bin/bash

#BSUB -o %J_STDOUT.log
#BSUB -e %J_STDERR.log
#BSUB -J Cstacks_${j}
#BSUB -n 3
#BSUB -M 2000000
#BSUB -q priority

## Simple script for submitting pstacks jobs

wd=/scratch/beegfs/monthly/cmatthey/data ## the directory containing the pstacks files
MM=1

mkdir -p $wd/data/cstacks/mm-$MM


samp=""
for i in $wd/pstacks/covmin-4/*allel*;
do
    sam+="-s ${i%%.*}* "
done;

module add UHTS/Analysis/stacks/1.30

cstacks -o $wd/cstacks/mm-$MM -n $MM -p 3 $sam
