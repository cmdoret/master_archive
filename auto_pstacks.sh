 

for M in 1 2 3 4 5 6;
do
    mkdir -p data/stacks/covmin-$M
    sed -i "s/\(M=\)[0-9]*/\1$M/g" src/stacks_pipeline/pstacks_script.sh
    bsub <./src/stacks_pipeline/pstacks_script.sh
done;
