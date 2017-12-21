# The purpose of this script is to generate a large table with one row per family, containing the following infos:
# FamID / Mother genotype at all CSD loci / Number of daughters / Proportion of het daugters at all CSD loci
# When mother genotypes are not available, infer them from offspring
# Cyril Matthey-Doret 
# 17.12.2017

#=== LOADING DATA ===#
pack <- c("ggplot2","dplyr", "readr", "viridis", "reshape2", "tibble")
lapply(pack, require, character.only = TRUE)
# Loading E/O/M genotype matrix
hom_path <- "data/assoc_mapping/grouped_keep_geno_EOM.tsv"
sum_stat <- read_tsv(file=hom_path, col_names=T)
# Loading list of individuals with family, generation and sex informations
indv <- read.table("data/individuals.tsv", header=T, sep='\t', stringsAsFactors = F)
# Loading list of SNPs with significant p-values for CSD association mapping
hits <- read.table("data/assoc_mapping/case_control/case_control_hits.tsv", header=T, stringsAsFactors = F)

#=== FORMATTING DATA ===#
# Filtering most significant hit on each chromosome
best <-  hits %>%
  group_by(Chr) %>%
  filter( fisher > 3 & fisher == max(fisher))

# Subset best CSD SNPs
csd_stat <- sum_stat %>%
  right_join(., best, by=c("Chr","BP")) %>%
  select(2,3,one_of(indv$Name[indv$Sex=='F']))


F3_geno <- csd_stat %>%
  select(1,2,one_of(indv$Name[indv$Generation=='F3']))

colnames(F3_geno[,colSums(F3_geno=='M')==0])
