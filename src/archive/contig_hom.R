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
colnames(contigs) <- c("Chr", "start","end")
# Loading 
hom_stat <- read_tsv(file='data/assoc_mapping/grouped_outpool_prophom.tsv', col_names=T, col_types = "icidddddd")

#==== PROCESS ===#

# Cleaning data file
hom_stat <- hom_stat[hom_stat$N.Males>0 & hom_stat$N.Females>0,]
hom_stat <- hom_stat[!is.na(hom_stat$Prop.Hom),]
hom_stat <- hom_stat[grep("chr.*",hom_stat$Chr),]
hom_stat$Chr <- gsub(x = hom_stat$Chr, pattern = "chr",replacement = "lf") 

# Computing proportion of homozygous individuals at each
  # 1. Group SNPs by contig (if overlap -> contig_row_number)
in_contig <- function(x){
  chr <- x[["Chr"]]; pos <- as.numeric(x[["BP"]])
  # Compares input position to set of genomic segments. Returns row index of the overlapping segment.
  # sub_tig <- contigs %>% filter(Chr==chr & start <= pos & end >= pos)
  tig_num <- with(contigs, which(as.character(Chr) == chr & pos >= start & pos <= end))
  return(tig_num)
}
hom_stat$tig <- apply(X = hom_stat, MARGIN = 1, FUN = in_contig)
  # 2. Compute mean homozygosity in each contig separately, weighting each SNP by sample size
hom_tig <- hom_stat %>% 
  group_by(tig) %>%
  summarise(weighted_hom=sum(N.Samples*Prop.Hom)/sum(N.Samples))
  # 3. Merging homozygosity values into contig table
out_contigs <- merge(hom_tig, contigs, by.x = "tig", by.y = "row.names")
#==== WRITE OUTPUT ===#

# Producing heatmap datafile for circos