 

for M in 1 2 3 4 5 6;
do
    bash ./src/stacks_pipeline/multi_pstacks.sh ./data/mapped/aln-4-19-100 $M
    for f in bsub_scripts/*;
    do
        bsub <./$f
    done;
done;
