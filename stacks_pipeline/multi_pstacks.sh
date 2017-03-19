#!/bin/bash

## This script will make a separate bsub script for each 
## sample in a directory. These scripts can then be 
## submitted with a for loop to the priority queue
## for super quick Ustacks analyses. See below for 
## arguments required:

wd=$1 ## the directory containing the fastq files
ID=$2 ## the stacks ID to start from


cd $wd


mkdir -p bsub_scripts
mkdir -p ustacks


for i in $(ls *fq.gz)
do 
    echo "Sample= $i, ID=$ID" 
    j=$(echo $i | cut -f1 -d '.')
    echo "#!/bin/bash" > ./bsub_scripts/bsub_${j}_script.sh
    echo "" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -L /bin/bash" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -o %J_STDOUT.log" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -e %J_STDERR.log" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -J Ustacks_${j}" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -n 3" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -M 2000000" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -q priority" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "" >> ./bsub_scripts/bsub_${j}_script.sh

    echo "module add UHTS/Analysis/stacks/1.30" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "" >> ./bsub_scripts/bsub_${j}_script.sh

    echo "ustacks -t gzfastq -f $(pwd)/$i -i $ID -r -d -p 3 -o $wd/ustacks" >> ./bsub_scripts/bsub_${j}_script.sh

    ((ID+=1))

done 
