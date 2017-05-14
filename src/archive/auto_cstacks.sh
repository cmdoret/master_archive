#!/bin/bash 

for M in 1 2 3 4;
do
    rm -fr data/cstacks/mm-$M;
    mkdir -p data/cstacks/mm-$M;
    sed -i "s/\(MM=\)[0-9]*/\1$M/g" src/stacks_pipeline/sub_cstacks.sh;
    bsub <./src/stacks_pipeline/sub_cstacks.sh;
done;
