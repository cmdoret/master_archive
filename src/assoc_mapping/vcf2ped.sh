# This scripts converts the VCF file containing all individuals to ped format.
# Takes a vcf file as input.
# Cyril Matthey-Doret
# 08.06.2017


vcftools --vcf $1 --plink  --out data/assoc_mapping/raw
