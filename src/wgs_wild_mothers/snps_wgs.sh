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


# Loading softwares
source src/misc/dependencies.sh
source src/mist/jobs_manager.sh

# Storing name and length of each contig in a file
cat $REF | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; }
                $0 !~ ">" {c+=length($0);}
                END { print c; }' > tig_sizes.tmp
# Removing empty lines
sed '/^\s*$/d' tig_sizes.tmp > ts.bak && mv ts.bak tig_sizes.tmp

# Progress tracker
tot_nucl=$(awk '{ sum += $2 } END { print sum }' tig_sizes.tmp)
done_nucl=0

# SNP calling is performed on regions  of 100kb in parallel to speed up the operation
# Iterating over contig
while read -a line
do

  tig_ID=${line[0]}
  tig_len=${line[1]}
  start=1;end=1
  # Splitting current contig into regions
  while [ $end -le $tig_len ]
  do

    # End takes smallest value between start+100kb and end of contig
    end=$(($tig_len<$end+100000?$tig_len:$end))
    region=$tig_ID":"$start"-"$end

    # Track progress
    let "done_nucl += $end-$start"
    prettyload $done_nucl $tot_nucl
# Do not submit more than 200 jobs at once
bmonitor WGSSNP 200
# 1 job per region
bsub << VAR
#!/bin/bash
#BSUB -L /bin/bash
#BSUB -o data/logs/wgs_snp-OUT.txt
#BSUB -e data/logs/wgs_snp-ERROR.txt
#BSUB -u cmatthey@unil.ch
#BSUB -J WGSSNP-$region
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q priority
#BSUB -R "rusage[mem=4000]"
#BSUB -M 4000000

# SNP calling in given region for all samples
samtools mpileup -f "$REF" \
                 -b <(find  "${WGS}/mapped/" -name "*sorted.bam" -type f ) \
                 -r "$region" \
                 -v > "${snps}/$region.tmp.vcf"
VAR

    let "start += 100000"

done < tig_sizes.tmp

rm tig_sizes.tmp

# Wait for all SNP calling jobs to finish before resuming
bmonitor WGSSNP 0

# Concatenating all regions VCF files into a large one
vcf-concat ${snps}/*.vcf > ${snps}/wild_mothers.vcf
# Removing all regions VCF
rm ${snps}/*.tmp.vcf
# Computing allelic diversity along genome
vcftools --window-pi $WIN --vcf ${snps}/wild_mothers.vcf --out nucleo_div
