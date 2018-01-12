#!/bin/bash
#BSUB -L /bin/bash
#BSUB -o data/logs/wgs_snp-OUT.txt
#BSUB -e data/logs/wgs_snp-ERROR.txt
#BSUB -u cmatthey@unil.ch
#BSUB -J WGSSNP
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -q normal
#BSUB -R "rusage[mem=64000]"
#BSUB -M 64000000

# Genome-wide variant calling and calculation of allelic diversity in a sliding window.
# Cyril Matthey-Doret
# 02.01.2018


# Parsing CL arguments
while [[ "$#" > 1 ]]; 
do 
    case $1 in
        # Working directory. Must contain folder with bam files
        --workdir) WGS="$2";;
        # Window size in which nucleotidic diversity is computed
        --winsize) WIN="$2";;
        # Reference genome
        --ref) REF="$2";;
        *) break;;
esac; shift; shift
done

module add UHTS/Analysis/vcftools/0.1.14
module add UHTS/Analysis/freebayes/0.9.9.2

snps="$WGS/variant/"
mkdir -p "$snps"
rm -rf "${snps}/*"

freebayes -f "$REF" "${WGS}/mapped/bam/*-uniq.sorted.bam" > "${snps}/wild_mothers.vcf"
vcftools --window-pi --vcf $WIN ${snps}/wild_mothers.vcf --out nucleo_div
