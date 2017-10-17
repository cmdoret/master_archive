# This script runs different steps to produce a linkage maps from processed reads.
# It follows the protocol described in Gadau, 2009 (doi:10.1101/pdb.prot5251)
# It is based on code snippets from Jens Bast.
# Cyril Matthey-Doret
# 16.10.2017


# Note that the reads must already be cleaned of adapters and trimmed.

# 1: Get clean 012 matrix

# Renaming fastq files for compatibility with MSTmap
declare -i int_val=1
for file in *.fq.gz;
do
  # Pad integer with 0s to have 3 digits
  printf -v pad_int "%03d" $int_val
  mv $file $pad_int".fq.gz"
  ((int_val++))
done

# Bowtie2 with RG flag included
for f in *.fq.gz
do
  bowtie2 --end-to-end \
          --very-sensitive \
          --rg-id $f \
          --threads 4 \
          -x ../mapping/Lys \
          -U $f \
          -S $f.sam
done

# Renaming RG flag to exclude file extensions from sample name
for align in *.sam;
do
  # Excluding file extension from name
  sample=${align%%.*}
  sed "s/${sample}.fq.gz/$sample/g" $align | \
  sed "s/\(@RG\tID:\)${sample}/\1${sample}\tSM:${sample}\tLB:1/g" > $align
done

# Convert to bam
for f in *.sam
do
  samtools view -bS -@ 24 $f -o $f.bam
done
# Sort bam files
for f in *.sam.bam
do
  samtools sort -@ 24 $f $f.sorted
done

# List of ploidy per sample: sample_name\t(1|2)
for file *.sorted.bam
do
  # Only haploid samples: set them all to 1 in the ploidy file
  echo $file | cut -d "." -f 1 | awk '{print$1,"\t1"}' | \
  sed 's/ //g' > names_ploidy
done

# Designed for samtools 0.1.19
samtools mpileup -u -D -q 20 \
                 -f ../canu2_low_corrRE.contigs.fasta \
                 -b files.txt | \
  bcftools view -c -g -v -s names_ploidy - > var.vcf

# Jens has vcftools 0.1.13
vcftools --vcf test_new.vcf \
         --minGQ 20 \
         --max-missing 0.8 \
         --maf 0.15 \
         --thin \
         --recode \
         --recode-INFO-all \
         --012 \
         --out test_new_filtered

# 2: Prepare MSTmap

# Requires an 'i' in front of sample names
cat Lys_var_filtered.012 | sed 's/^/i/'

# Transpose rows and col (find oneliner substitute to Jens' script)
~/Dropbox/scripts/transpose.sh test | \
  sed 's/ /\t/g' > Lys_var_filtered_012.tp

# Prepend loci names
cat Lys_var_filtered.012.pos | \
  sed 's/\t/_/g' | \
  sed '1s/^/locus_name\n/' > loci
paste loci Lys_var_filtered_012.tp > Lys_MSTmap_in

# Change encoding from 012 to UAB
cat Lys_MSTmap_in | \
  sed 's/\t-1/\tU/g' | \
  sed 's/\t0/\tA/g' | \
  sed 's/\t2/\tB/g' > Lys_MSTmap.A

# Double and reverse phase
cat Lys_MSTmap_in | \
  sed 's/\t-1/\tU/g' | \
  sed 's/\t0/\tB/g' | \
  sed 's/\t2/\tA/g' > Lys_MSTmap.B

# Reversed phase file and delete first ID line
# Locus ID will be A/B switched in name
cat Lys_MSTmap.B | sed '1d' > Lys_MSTmap.B2

# Test: Apply tag "RP_" to locus ID:
cat Lys_MSTmap.B | \
  sed '1d' | \
  sed 's/^/RP_/' > Lys_MSTmap.B3

# Join both
cat Lys_MSTmap.A Lys_MSTmap.B2 > tmp

# Add header after setting parameters in a txt file
cat input tmp > Lys_MSTmap.inA

# input parameters in extra file (input):

population_type DH
population_name Lys
distance_function kosambi
cut_off_p_value 1e-6
no_map_dist 15.0
no_map_size 2 #set to 0 to disable this
missing_threshold 1.00
estimation_before_clustering no
detect_bad_data yes
objective_function COUNT
number_of_loci 2209 #(double number of loci because doubled, .i.e: wc -l tmp (-1 header line))
number_of_individual 124 #(number of samples)

# 3: Running MSTmap
~/Software/MSTMap/MSTMap.exe Lys_MSTmap.inA Lys_MSTmap_p1e-6.out
# cut_off_p_value 1e-6
# gives 6 linkage groups
