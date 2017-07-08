#!/bin/bash

## This script will make a separate bsub script for each 
## sample in a directory. These scripts can then be 
## submitted with a for loop to the priority queue
## for super quick pstacks analyses. See below for 
## arguments required:

wd=$1 ## the directory containing the sam files
declare -i ID=1 ## the stacks ID to start from

rm -rf bsub_scripts
mkdir -p bsub_scripts
rm -rf data/pstacks/covmin-$2
mkdir -p data/pstacks/covmin-$2
rm -rf data/logs/pstacks
mkdir -p data/logs/pstacks


# Writing scripts automatically
for i in $(ls $wd/bam/*uniq.sorted.bam)
do
    echo "Sample= $i, ID=$ID" 
    j=$(echo ${i##*/} | cut -f1 -d '.')
    echo "#!/bin/bash" > ./bsub_scripts/bsub_${j}_script.sh
    echo "" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -L /bin/bash" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -o ./data/logs/pstacks/PST_COVMIN$2_${j}_STDOUT.log" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -e ./data/logs/pstacks/PST_COVMIN$2_${j}_STDERR.log" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -J PST${j}" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -n 3" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -M 2000000" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -q priority" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "" >> ./bsub_scripts/bsub_${j}_script.sh

    echo "module add UHTS/Analysis/stacks/1.30" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "" >> ./bsub_scripts/bsub_${j}_script.sh

    echo "pstacks -f $i -i $ID -o ./data/pstacks/covmin-$2 -m $2 -p 3 -t bam" >> ./bsub_scripts/bsub_${j}_script.sh

    ID+=1
 
done


# Submitting all scripts in a loop
for f in bsub_scripts/*;
do
    bsub <./$f
    sleep 1;
done;

# Waiting until all Pstacks jobs are finished
while [ $(bjobs -w | awk '/RUN/ {print $7}' | grep 'PST' | wc -l) -gt 0 ]
do
    sleep 2;
done

# Cleaning filenames to match popmap later on
for f in ./data/pstacks/covmin-$2/*.tsv*; 
do 
    mv $f $(echo "$f" | sed 's/-BWA-uniq.sorted//');
done;
