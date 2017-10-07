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

# Map families to mother categories
sum_stat$cluster <- groups$cluster[match(sum_stat$Family, groups$Family)]
sum_stat <- sum_stat[!is.na(sum_stat$cluster),]

# Group SNP by mother category
cat_stat <- sum_stat %>%
  group_by(Locus.ID, Chr, BP, Family) %>%
  summarise(Nf=sum(N.Females), Nm=sum(N.Males), N=sum(N.Samples), 
            Prop.Hom=sum(Prop.Hom*N.Samples, na.rm=T)/sum(N.Samples, na.rm=T), 
            Prop.Hom.F=sum(Prop.Hom.F*N.Females, na.rm=T)/sum(N.Females, na.rm=T), 
            Prop.Hom.M=sum(Prop.Hom.M*N.Males, na.rm=T)/sum(N.Males, na.rm=T))

# Compute assocation on each cat. separately
# Output single file with corrected significant p-values with
# associatied SNP and category.

#==== COMPUTE STATISTICS ====#
# Abbreviations: M: Males, F: Females, T: Male+Female, 
# o:homozygous, e:heterozygous, t:hom+het, E: Expected

get_fisher <- function(df){
  mat <- matrix(as.numeric(df[c(8:11)]), ncol=2)
  f <- fisher.test(as.table(mat), alt="two.sided")
  return(f$p.value)
}


odds_list <- cat_stat %>% 
  rename(Ft = Nf, Mt = Nm, Tt = N) %>%
  mutate(Fo = Ft * Prop.Hom.F, Mo = Mt * Prop.Hom.M, 
         Fe = Ft * (1-Prop.Hom.F), Me = Mt * (1-Prop.Hom.M),
         To = Tt * Prop.Hom, Te = Tt * (1-Prop.Hom)) %>%
  select(-Prop.Hom, -Prop.Hom.F, -Prop.Hom.M) %>%
  mutate_at(funs(round(.,0)), .vars = c("Fo","Fe","Mo","Me"))

odds_list$fisher <- apply(odds_list, 1,  get_fisher)

#odds_list$fisher <- p.adjust(odds_list$fisher, method = "bonferroni")
for(group in unique(odds_list$Family)){
  odds_list$fisher[odds_list$Famly==group] <- p.adjust(odds_list$fisher[odds_list$Family==group], method = "BH")
}
nloci <- log2(max(groups$cluster)+1)
#========= VISUALISE ========#
odds_chrom <- odds_list[grep("chr.*",odds_list$Chr),]
pdf(paste0(out_folder, "/../plots/","case_control_hits_",nloci,"loci.pdf"), width=12, height=12)
ggplot(data=odds_chrom, aes(x=BP, y=-log10(fisher), col=Family)) + geom_point() + facet_grid(~Chr, scales='free_x') +  
  geom_hline(aes(yintercept=-log10(0.05))) + geom_hline(aes(yintercept=-log10(0.01)), lty=2, col='red') + 
  xlab("Genomic position") + ylab("-log10 p-value") + ggtitle("Case-control associaiton test for CSD") + ylim(c(0,10))
dev.off()
#======= WRITE OUTPUT =======#
# Number of groups is (2^n)-1 where n is the number of CSD loci
nloci <- log2(max(groups$cluster)+1)
write.table(odds_chrom, paste0(out_folder, "case_control_hits_", nloci , "loci.tsv"), 
            sep='\t', row.names=F, quote=F)