# This script parses the populations output vcf file to analyse and visualise coverage across the genome.
# Cyril Matthey-Doret
# 27.07.2017

library(vcfR) 
library(ape)
# Simple VCF file import and parsing

vcf_fam <- read.vcfR("dummy.vcf")
ref_genome <- read.dna("dummy.fasta", format='fasta')
chrom_test <- create.chromR(vcf=vcf_fam,seq=ref_genome)


for(chrom in names(ref_genome)){

  

plot(chrom_test)

myDNA <- vcfR2DNAbin(vcfR_test, ref.seq = myRef)
seg.sites(myDNA)
image(myDNA)
