cd "$(dirname "$0")"
mkdir -p vcftools

for folder in ../../data/populations/*;
do
    outname=$(basename $folder);
    vcftools --vcf $folder/*.vcf --out ./vcftools/$outname --depth;
    vcftools --vcf $folder/*.vcf --out ./vcftools/$outname --het;
    paste ./vcftools/$outname.het ./vcftools/$outname.idepth \
     | cut -f6,7 --complement > vcftools/summary_$outname;
done;
