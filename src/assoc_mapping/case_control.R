# In this script I use a naive test for allelic association
# to test for association with homozygosity in males. The test should be performed 
# only within a category of family (inferred from diploid male production).
# Cyril Matthey-Doret
# 18.08.2017

#======== LOAD DATA =========#
library(dplyr);library(readr);library(ggplot2)
hom_path <- "../../data/assoc_mapping/grouped_outpool_prop_hom_fixed_sites.tsv"
sum_stat <- read_tsv(file=hom_path, col_names=T, col_types = "iciddddddc")

#======= PROCESS DATA =======#
sum_stat <- sum_stat[sum_stat$N.Males>0 & sum_stat$N.Females>0,]
sum_stat <- sum_stat[!is.na(sum_stat$Prop.Hom),]

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
