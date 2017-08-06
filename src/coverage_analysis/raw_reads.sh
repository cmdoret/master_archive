# This script stores the number of processed reads per individuals in a table.
# Cyril Matthey-Doret
# 06.08.2017

read_dir="../../data/processed/"
indv="../../data/individuals"
out="../../data/coverage/reads_per_sample.txt"
echo "sample,reads" > $out

for sample in $(cut -f1 $indv | tail -n +2 | sort | uniq);
do
    count=$(zgrep "@" $read_dir/$sample.fq.gz | wc -l)
    echo $sample","$count >> $out
done
