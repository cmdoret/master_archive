
# This script is used for copying all snps, tags and alleles files from valid pstacks samples to sstacks folder
# It takes 3 arguments: the paths of P, C and S stacks folders, respectively
PSTACK=$1
CSTACK=$2
SSTACK=$3
for f in $SSTACK/*; 
do
    tf=$(basename ${f%%.matches*}); 
    cp $PSTACK/$tf* $SSTACK; \
done;
cp $CSTACK/* $SSTACK;