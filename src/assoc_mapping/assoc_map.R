# In this script I use a modified version of the basic test for allelic association
# to test for association with homozygosity in males. The test should be performed 
# only within a category of family (inferred from diploid male production).
# Cyril Matthey-Doret
# 18.08.2017

#======== LOAD DATA =========#
library(dplyr);library(readr)
hom_path <- "../../data/assoc_mapping/grouped_prop_hom_fixed_sites.tsv"
sum_stat <- read_tsv(file=hom_path, col_names=T, col_types = "iciddddddc")
sum_stat <- sum_stat[sum_stat$N.Males>0 & sum_stat$N.Females>0,]
sum_stat <- sum_stat[!is.na(sum_stat$Prop.Hom),]

#======= PROCESS DATA =======#


#==== COMPUTE STATISTICS ====#


#========= VISUALISE ========#


#======= WRITE OUTPUT =======#
