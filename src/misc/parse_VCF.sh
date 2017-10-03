# Parsing the output VCF file from populations runs in each family
# into summarized tables and convenient genotype matrices.
# Takes 2 argument: input path and whether populations was run separately 
# for each family (F), or once for all individuals (T)
# Cyril Matthey-Doret
# 24.06.2017

if [[ $# -eq 0 ]]; then
    echo "No argument given, exiting now."
    echo "I need the folder where populations output files are stored."
    exit 0
fi
declare -i increm=0
out_path=data/ploidy/vcftools
rm -rf $out_path
mkdir -p $out_path

group=$2  
# did populations run separately for each family (F) or only once for all individuals (T).

if [ $group = "T" ]
then
    mkdir -p $out_path/all/
    if [ -f $1/*.vcf ]
    then        
        vcftools --vcf $1/*.vcf --out $out_path/all/all --depth;
        # Mean read depth per individual
        vcftools --vcf $1/*.vcf --out $out_path/all/all --het;
        # Mean heterozygosity per individual
        vcftools --vcf $1/*.vcf --out $out_path/all/all --012;
        # Genotype matrix[SNP, individual]
        paste $out_path/all/all.het $out_path/all/all.idepth \
        | cut -f-5,8- > $out_path/summary_full.txt;
        # Putting together mean depth and mean heterozygosity
    fi
else
    for folder in $1/*;
    # Iterating over families
    do
        if [ -f $folder/*.vcf ];
        then
            fam=$(basename $folder);  # Name of family
            mkdir -p $out_path/$fam  # Creating one separate folder per family
            vcftools --vcf $folder/*.vcf --out $out_path/$fam/$fam --depth;
            # Mean read depth per individual
            vcftools --vcf $folder/*.vcf --out $out_path/$fam/$fam --het;
            # Mean heterozygosity per individual
            vcftools --vcf $folder/*.vcf --out $out_path/$fam/$fam --012;
            # Genotype matrix[SNP, individual]
            paste $out_path/$fam/$fam.het $out_path/$fam/$fam.idepth \
            | cut -f-5,8- > $out_path/$fam/summary_$fam;
            # Putting together mean depth and mean heterozygosity
            if [[ $increm -eq 0 ]]; then
                ((increm++))
                head -1 $out_path/$fam/summary_$fam > $out_path/summary_full.txt
                # Creating summary file common to all families and adding header
            fi
            tail -n +2 $out_path/$fam/summary_$fam >> $out_path/summary_full.txt
            # Appending each family-specific het + depth summary to the common one
            # without the header
        fi
    done
fi
sed -i'' '/nan/d' $out_path/summary_full.txt
