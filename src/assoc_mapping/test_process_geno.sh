#BSUB -q normal
#BSUB -M 20000000
#BSUB -R "rusage[mem=20000]"
#BSUB -n 8
#BSUB -R "span[ptile=8]"
group="T"
pop="data/populations/grouped_d-3_r-80/"
out="data/assoc_mapping/"
python2 process_genomic.py $pop $out $group
