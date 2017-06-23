
# This script is used for copying all snps, tags and alleles files from valid pstacks samples to sstacks folder
# It takes 3 arguments: the paths of P, C and S stacks folders, respectively

PSTACK=$1
CSTACK=$2
SSTACK=$3

if [[ $# -eq 0 ]] ; then
    echo 'I need 3 arguments ! The path to pstacks, cstacks and sstacks output files, respectively. Exiting.'
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
for fam in $(cut -f3 data/individuals|sort|tail -n +2|uniq)
# All families in dataset (excluding header with tail)
do
    mkdir -p $SSTACK/$fam
    for indv in $(awk -v var="$fam" '$3 == var {print $1}' data/individuals)
    # All individuals in each family
    do
        mv $SSTACK/$indv* $SSTACK/$fam  # Moving all files of each individual to its family folder
        cp $SSTACK/*catalog* $SSTACK/$fam/  # Copying the whole catalog for each family
    done
    echo "Family $fam folder ready for populations."
done    
