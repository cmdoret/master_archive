 
# This script generates files containing average coverage per site in each family.
# The files are then merged into 1 table containing all families.
# Cyril Matthey-Doret
# 27.07.2017

indv='../../data/individuals'  
# Contains various information about samples
pop='../../data/populations/d-3_r-80'
# Contains output from stacks populations module
out='../../data/coverage/'
# Desired output folder
families=$(cut -f3 $indv | tail -n +2 | sort | uniq)
# Extracting family info
mkdir -p $out  # creating output directory
grouped="T"  # Grouped populations run (T) or 1 per family (F)


if [ $grouped = "F" ]  # If one populations run per family
then
    for fam in $families;  # iterate over families
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
else  # One pooled populations run
    vcftools --site-mean-depth --vcf $pop/*.vcf --out $out/all
    rm $out/all.log
    # Generating read depth per site averaged across samples and 
    # removing log files
    echo -e $(head -n 1 $out/all.ldepth.mean)"\tFamily" > $out/site_depth.txt
    # Creating new output file and including header
    tail -n +2 $out/all.ldepth.mean > tmp1_all
    # extracting each family's file content without header
    awk 'BEGIN{OFS="\t"}{print $0, "all";}' tmp1_all > tmp2_all
    # Adding column with family information 
fi

cat tmp2_* >> $out/site_depth.txt
rm tmp*