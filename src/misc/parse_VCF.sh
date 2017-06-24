# Parsing the output VCF file from populations runs in each family
# into summarized tables and convenient genotype matrices.
# Cyril Matthey-Doret
# 24.06.2017

if [[ $# -eq 0 ]]; then
    echo "No argument given, exiting now."
    echo "I need the folder where populations output files are stored."
    exit 0
fi

out_path=data/ploidy/vcftools
mkdir -p $out_path

for folder in $1/*;
do
    outname=$(basename $folder);
    mkdir -p $out_path/$outname
    vcftools --vcf $folder/*.vcf --out $out_path/$outname/$outname --depth;
    vcftools --vcf $folder/*.vcf --out $out_path/$outname/$outname --het;
    vcftools --vcf $folder/*.vcf --out $out_path/$outname/$outname --012;
    paste $out_path/$outname/$outname.het $out_path/$outname/$outname.idepth \
     | cut -f6,7 --complement > $out_path/$outname/summary_$outname;
done;
