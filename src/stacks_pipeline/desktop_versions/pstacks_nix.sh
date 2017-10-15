#!/bin/bash

# This script will run pstacks in a loop without using an LSF cluster.
# Cyril Matthey-Doret
# 10.10.2017


declare -i ID=1 ## the stacks ID to start from


# parsing CL arguments
while [[ "$#" > 1 ]]; do case $1 in
    # Mimimum stack coverage for pstacks
    --m) M="$2";;
    # Path to the folder containing input alignment files
    --map) map="$2";;
    # Path to output stacks files
    --out) out_dir="$2";;
    # Path to log folder
    --log) logs="$2";;
    *) break;;
  esac; shift; shift
done

# Cleaning directories
rm -rf $out_dir; mkdir -p $out_dir
rm -rf $logs/pstacks; mkdir -p $logs/pstacks

# Writing scripts automatically
for i in $(ls $wd/bam/*uniq.sorted.bam)
do
    echo "Sample= $i, ID=$ID"
    j=$(echo ${i##*/} | cut -f1 -d '.')
    pstacks -f $i -i $ID -o $out_dir -m $M -p 3 -t bam

    ID+=1

done

# Cleaning filenames to match popmap later on
for f in $out_dir/*.tsv*;
do
    mv $f $(echo "$f" | sed 's/-BWA-uniq.sorted//');
done;
