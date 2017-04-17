#!/bin/bash

## This script will make a separate bsub script for each 
## sample in a directory. These scripts can then be 
## submitted with a for loop to the priority queue
## for super quick cstacks analyses. See below for 
## arguments required:

wd=$1 ## the directory containing the pstacks files


mkdir -p bsub_scripts
mkdir -p data/cstacks/mm-$2


for i in $(ls $wd/CF*01*.allele*)
do
    echo "Sample= $i" 
    j=$(echo ${i##*/} | cut -f1 -d '_')
    echo "#!/bin/bash" > ./bsub_scripts/bsub_${j}_script.sh
    echo "" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -L /bin/bash" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -o %J_STDOUT.log" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -e %J_STDERR.log" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -J Cstacks_${j}" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -n 3" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -M 2000000" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -q priority" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "" >> ./bsub_scripts/bsub_${j}_script.sh

    echo "module add UHTS/Analysis/stacks/1.30" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "" >> ./bsub_scripts/bsub_${j}_script.sh

    echo "cstacks -s ${i%%.*} -o ./data/cstacks/mm-$2 -n $2 -p 3" >> ./bsub_scripts/bsub_${j}_script.sh

done
