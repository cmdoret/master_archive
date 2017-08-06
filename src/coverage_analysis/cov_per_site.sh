 
# This script generates files containing average coverage per site in each family.
# The files are then merged into 1 table containing all families.
# Cyril Matthey-Doret
# 27.07.2017

indv='../../data/individuals'  
# Contains various information about samples
pop='../../data/populations/d-5_r-10'
# Contains output from stacks populations module
out='../../data/coverage/'
# Desired output folder
families=$(cut -f3 $indv | tail -n +2 | sort | uniq)
# Extracting family info
mkdir -p $out  # creating output directory

for fam in $families;
do
    vcftools --site-mean-depth --vcf $pop/$fam/*.vcf --out $out/$fam
    rm $out/$fam.log
    # Generating read depth per site averaged across samples and 
    # removing log files
done

add_head=0
for fam in $families;
do
    if [ $add_head -eq 0 ];  # Only runs on first file
    then
        echo -e $(head -n 1 $out/$fam.ldepth.mean)"\tFamily" > $out/site_depth.txt
        # Creating new output file and including header
        add_head=1
    fi
    tail -n +2 $out/$fam.ldepth.mean > tmp1_$fam
    # extracting each family's file content without header
    awk -v var=$fam 'BEGIN{OFS="\t"}{print $0, var;}' tmp1_$fam > tmp2_$fam
    # Adding column with family information 
done

cat tmp2_* >> $out/site_depth.txt
rm tmp*
