#!/bin/bash
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

snps="$WGS/variant/"
mkdir -p "$snps"
rm -rf "${snps}/*"

bsub << VAR
#!/bin/bash
#BSUB -L /bin/bash
#BSUB -o data/logs/wgs_snp-OUT.txt
#BSUB -e data/logs/wgs_snp-ERROR.txt
#BSUB -u cmatthey@unil.ch
#BSUB -J WGSSNP
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -q long
#BSUB -R "rusage[mem=16000]"
#BSUB -M 16000000

# Loading softwares
source src/misc/dependencies.sh

#find  "${WGS}/mapped/" -name "*sorted.bam" -type f | \
#    xargs freebayes -f "$REF" > "${snps}/wild_mothers.vcf"

samtools mpileup -f "$REF" \
                 -b <(find  "${WGS}/mapped/" -name "*sorted.bam" -type f ) \
                 -v > "${snps}/wild_mothers.vcf"

vcftools --window-pi $WIN --vcf ${snps}/wild_mothers.vcf --out nucleo_div
VAR
