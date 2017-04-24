 
#!/bin/bash

## This script runs sstacks, matching samples against the catalogue.
## It will exclude samples with very low number of reads from the analysis (less than 10% of
## the mean number of reads across samples), as they were not used when building the catalogue

wd=/scratch/beegfs/monthly/cmatthey/data ## the directory containing the pstacks files
declare -i n=0 tot=0


rm -rf bsub_scripts/
mkdir -p bsub_scripts
rm -rf data/sstacks
mkdir -p data/sstacks

for f in $wd/pstacks/covmin-4/*tags*;
do
    tot+=$(zcat $f | wc -l);
    n+=1;
done;

samp=""
for i in $wd/pstacks/covmin-4/*tags*;
do
    if [ "$(zcat $i | wc -l)" -gt $(($tot/($n*10))) ];
    then
        samp+="${i%%.tags*} "
    fi;
done;

for i in $samp
do
    echo "Sample= $i, ID=$ID" 
    j=$(echo ${i##*/} | cut -f1 -d '.')
    echo "#!/bin/bash" > ./bsub_scripts/bsub_${j}_script.sh
    echo "" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -L /bin/bash" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -o %J_STDOUT.log" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -e %J_STDERR.log" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -J SST${j}" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -n 3" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -M 2000000" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "#BSUB -q priority" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "" >> ./bsub_scripts/bsub_${j}_script.sh

    echo "module add UHTS/Analysis/stacks/1.30" >> ./bsub_scripts/bsub_${j}_script.sh
    echo "" >> ./bsub_scripts/bsub_${j}_script.sh

    echo "sstacks -b 1 -c $1/batch_0 -s $i -o ./data/sstacks/ -p 3" >> ./bsub_scripts/bsub_${j}_script.sh

done

for f in bsub_scripts/*;
do
    bsub <./$f;
done;

while [ $(bjobs -w | awk '/RUN/ {print $7}' | grep 'SST' | wc -l) -gt 0 ]
do
    sleep 2;
done;
