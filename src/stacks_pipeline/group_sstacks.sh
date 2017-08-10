
# This script is used for copying all snps, tags and alleles files from valid pstacks samples to sstacks folder
# It takes 4 arguments: the paths of P, C and S stacks folders, respectively, and whether (T) or not (F) it should 
# separate STACKS files into subfolders by family.
# Cyril Matthey-Doret
# 08.08.2017

PSTACK=$1
CSTACK=$2
SSTACK=$3
group=$4

if [[ $# -eq 0 ]] ; then
    echo 'I need 4 arguments ! The path to pstacks, cstacks and sstacks output files, respectivels, and a boolean (T/F) to know whether or not all samples will be grouped in a populations run. Exiting.'
    exit 0
fi

rm -f $SSTACK/batch*  # Remove files from previous runs
for f in $SSTACK/*;   # Iterating over sstacks files
do
    tf=$(basename ${f%%.matches*});   # Extracting individual names
    cp $PSTACK/$tf* $SSTACK   # Copying pstacks files corresponding to names to sstacks folder
done;
cp $CSTACK/* $SSTACK;  # Copying catalogue as well

# Below: organizing subfolders to run populations separately on each family
if [ $group = "F" ]
then

    for fam in $(cut -f3 data/individuals|tail -n +2|sort|uniq)
    # All families in dataset (excluding header with tail)
    do
        mkdir -p $SSTACK/$fam
        for indv in $(awk -v var="$fam" '$3 == var {print $1}' data/individuals)
        # All individuals in each family
        do
            mv $SSTACK/$indv* $SSTACK/$fam 2> /dev/null || echo "Individual $indv excluded from the analysis."
            # Moving all files of each individual to its family folder
        done
        cp $SSTACK/*catalog* $SSTACK/$fam/  # Copying the whole catalog for each family
        echo "Family $fam folder ready for populations."
    done    
fi