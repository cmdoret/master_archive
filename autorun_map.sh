
MM=4
ALG=mem
for K in 15 17 19 21;
do
    for W in 90 100 110 120;
    do
        mkdir -p data/mapping/$ALG-$MM-$K-$W
        sed -i "s/\(MM=\)[0-9]*/\1$MM/g" src/mapping/bwa_script.sh
        sed -i "s/\(ALG=\)[a-z]*/\1$ALG/g" src/mapping/bwa_script.sh
        sed -i "s/\(K=\)[0-9]*/\1$K/g" src/mapping/bwa_script.sh
        sed -i "s/\(W=\)[0-9]*/\1$W/g" src/mapping/bwa_script.sh
        bsub <./src/mapping/bwa_script.sh
    done;
done;
K=19
W=100
ALG=aln
for MM in 0 2 4 6 8;
do
    mkdir -p data/mapping/$ALG-$MM-$K-$W
    sed -i "s/\(MM=\)[0-9]*/\1$MM/g" src/mapping/bwa_script.sh
    sed -i "s/\(ALG=\)[a-z]*/\1$ALG/g" src/mapping/bwa_script.sh
    sed -i "s/\(K=\)[0-9]*/\1$K/g" src/mapping/bwa_script.sh
    sed -i "s/\(W=\)[0-9]*/\1$W/g" src/mapping/bwa_script.sh
    bsub <./src/mapping/bwa_script.sh
done;
