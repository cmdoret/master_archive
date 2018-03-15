# Computes PI nucleotidic diversity from WGS samples at all SNPs from the RADseq analysis.
# Cyril Matthey-Doret
# 14.03.2018

source src/misc/dependencies.sh

# Get list of sites from association mapping files
sites="data/assoc_mapping/case_control/case_control_all.tsv"
wgs_vcf="data/wgs_wild/variant/wild.sorted.vcf"
out_dir="data/wgs_wild/RAD_sites"
mkdir -p $out_dir

# Generate site list for vcftoolsa
tail -n +2 $sites | awk 'BEGIN{OFS="\t"}{print $2,$3}' > "$out_dir/sites_list.tsv"

# Compute PI per site with vcftools
vcftools --positions $out_dir/sites_list.tsv --site-pi --vcf $wgs_vcf --out $out_dir/rad2wgs_pi.vcf
