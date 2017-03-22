#!bin/bash 

MAL=$(awk -F, '/M/ {print $1}' F4_table.csv | cat)
FEM=$(awk -F, '/F/ {print $1}' F4_table.csv | cat)

files=*.fq
for m in $MAL;
do
    for i in $files;
    do
        if [ $i == $m'.fq' ];
        then
            mv $i $m'_M.fq'
        fi
    done
        
done

for f in $FEM;
do
    for i in $files;
    do
        if [ $i == $f'.fq' ];
        then
            mv $i $f'_F.fq'
        fi
    done
done
