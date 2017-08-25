# In this script I use a naive test for allelic association
# to test for association with homozygosity in males. The test should be performed 
# only within a category of family (inferred from diploid male production).
# Cyril Matthey-Doret
# 18.08.2017

#======== LOAD DATA =========#
library(dplyr);library(readr);library(ggplot2)
in_args <- commandArgs(trailingOnly = T)
groups_path <- in_args[1]
groups <- read.table(groups_path, header = T)
hom_path <-in_args[2]
sum_stat <- read_tsv(file=hom_path, col_names=T, col_types = "iciddddddc")
out_folder <- in_args[3]

#======= PROCESS DATA =======#
sum_stat <- sum_stat[sum_stat$N.Males>0 & sum_stat$N.Females>0,]
sum_stat <- sum_stat[!is.na(sum_stat$Prop.Hom),]

#!!!!!! TODO !!!!!!#
# Map families to mother categories
# Group SNP by mother category
# Compute assocation on each cat. separately
# Output single file with corrected significant p-values with
# associatied SNP and category.

#==== COMPUTE STATISTICS ====#
# Abbreviations: M: Males, F: Females, T: Male+Female, 
# o:homozygous, e:heterozygous, t:hom+het, E: Expected
odds_list <- sum_stat %>% 
  rename(Ft = N.Females, Mt = N.Males, Tt = N.Samples) %>%
  mutate(Fo = Ft * Prop.Hom.F, Mo = Mt * Prop.Hom.M, 
         Fe = Ft * (1-Prop.Hom.F), Me = Mt * (1-Prop.Hom.M),
         To = Tt * Prop.Hom, Te = Tt * (1-Prop.Hom)) %>%
  select(-Prop.Hom, -Prop.Hom.F, -Prop.Hom.M) %>% 
  mutate(EFo = Ft * To, EMo = Mt * To, EFe = Ft * Te, EMe = Mt * Te) %>%
  mutate(chi2 = (Fo-EFo)^2/EFo + (Mo-EMo)^2/EMo + (Fe-EFe)^2/EFe + (Me-EMe)^2/EMe)

odds_list$pval <- sapply(X = odds_list$chi2, function(x) pchisq(x, df=1, lower.tail=F, log.p=T))

#========= VISUALISE ========#
odds_chrom <- odds_list[grep("chr.*",odds_list$Chr),]
ggplot(data=odds_chrom, aes(x=BP, y=abs(pval))) + facet_grid(~Chr, scales='free_x') + geom_point() + 
  xlab("Genomic position") + ylab("-log2 p-value") + ggtitle("Case-control associaiton test for CSD")

#======= WRITE OUTPUT =======#
# Number of groups is (2^n)-1 where n is the number of CSD loci
nloci <- log2(max(groups$cluster)+1)
write.table(odds_chrom, paste0(out_path, "case_control_hits_", nloci , "loci.tsv"), 
            sep='\t', row.names=F, quote=F)