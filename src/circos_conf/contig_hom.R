# This script is used solely to generate a heatmap data file 
# for circos. It computes percentage of homozygosity relative 
# to mother, a proxy of recombination rate in central-fusion 
# automixis, in each contig separately.
# Cyril Matthey-Doret
# 23.11.2017

#==== SETUP ENVIRONMENT ===#

# Loading libraries for faster data loading and polyvalent manipulations
library(readr);library(dplyr)
setwd("~/Documents/UNIL/Master/master_project/")

# Loading coordinates of original contigs anchored into chromosomes
contigs <- read.table(file='data/circos/tig_lim.lf.txt')
# Loading 
hom_stat <- read_tsv(file='data/assoc_mapping/grouped_outpool_prophom.tsv', col_names=T, col_types = "icidddddd")

#==== PROCESS ===#

# Cleaning data file
hom_stat <- hom_stat[hom_stat$N.Males>0 & hom_stat$N.Females>0,]
hom_stat <- hom_stat[!is.na(hom_stat$Prop.Hom),]

# Computing proportion of homozygous individuals at each
  # 1. Group SNPs by contig (if overlap -> contig_row_number)
  # 2. Compute mean homozygosity in each contig separately, weighting each SNP by sample size

#==== WRITE OUTPUT ===#

# Producing heatmap datafile for circos