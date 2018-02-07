#!/bin/bash

# This script runs one ustacks job per sample in a loop on an LSF
# cluster. This speeds up the process by heavily parallelizing tasks.
# Cyril Matthey-Doret
# 10.10.2017

source src/misc/jobs_manager.sh
declare -i ID=1 ## the stacks ID to start from

# parsing CL arguments
while [[ "$#" > 1 ]]; do case $1 in
    # Mimimum stack coverage for ustacks
    --m) M="$2";;
    # Maximum nucleotide distance allowed between stacks
    --mm) MM="$2";;
    # Path to the folder containing input processed reads
    --reads) reads="$2";;
    # Path to output stacks files
    --out) out_dir="$2";;
    # Path to log folder
    --log) logs="$2";;
    *) break;;
  esac; shift; shift
done

# Cleaning directories
rm -rf $out_dir; mkdir -p $out_dir
rm -rf $logs/ustacks; mkdir -p $logs/ustacks

# For each sample
for i in $(ls $reads/*.fq*)
do
    # Do not queue more than 100 jobs at a time
    bmonitor UST 100
    echo "Sample= $i, ID=$ID"
    j=$(echo ${i##*/} | cut -f1 -d '.')

    bsub <<UST
    #!/bin/bash
    #BSUB -L /bin/bash
    #BSUB -o $logs/ustacks/UST_COV${M}MM${MM}_${j}_STDOUT.log
    #BSUB -e $logs/ustacks/UST_COV${M}MM${MM}_${j}_STDERR.log
    #BSUB -J UST${j}
    #BSUB -n 3
    #BSUB -M 2000000
    #BSUB -q priority

    # Loading softwares
    source src/misc/dependencies.sh
    ustacks -f $i -i $ID -o $out_dir -m $M -M $MM -p 3
UST
    ID+=1
done

# Waiting until all Pstacks jobs are finished
bmonitor UST 0
