# This script is used to differentiate between telocentric and metacentric
# chromosomes. This is done by using the density of heterozygous sites in the
# mother that become homozygous in offspring along the chromosome.
# 30.07.2017
# Cyril Matthey-Doret

library(ggplot2); library(dplyr)


# Genome statistics
stat_path <- '../../data/populations/d-20_r-80/'
indv <- read.table('../../data/individuals', header=T)
# phi_path <- commandArgs(TrailingOnly=T)[1]

# removing SNPs homozygous in all females (consequently in mothers)
sum_stat <- data.frame()
for(fam in list.dirs(stat_path)[2:length(list.dirs(stat_path))]){  # Excluding first dir (parent)
  tmp_stat <- read.csv(paste0(fam,'/batch_0.sumstats.tsv'),header=T,skip=2,sep='\t')
  tmp_stat$fam <- rep(basename(fam))
  sum_stat <- rbind(sum_stat, tmp_stat)
  mother_hom <- NULL
}
rm_snps<- sum_stat %>% 
  group_by(fam,Locus.ID,Col) %>%  # Each group contains one SNP (both male and female pop)
  summarise(hom_mot = sum(Pop.ID=='F' & Q.Nuc=='-'))  # If all females are homozygous -> mother is homozygous

rm_snps <- rm_snps[rm_snps$hom_mot==0, 1:3]
sum_stat <- merge(sum_stat,rm_snps,by=c("Locus.ID","Col","fam"))

# Computing total homozygousity and number of individuals per SNP
male_stat <- sum_stat[sum_stat$Pop.ID=="M",c('Chr','BP','Obs.Hom', 'fam', 'N',"Locus.ID")]
female_stat <- sum_stat[sum_stat$Pop.ID=="F",c('Chr','BP','Obs.Hom', 'fam', 'N')]
chrom_stat <- merge(male_stat,female_stat,by=c('Chr','BP','fam'))
attach(chrom_stat)
chrom_stat$hom <- (Obs.Hom.x * N.x + Obs.Hom.y * N.y)/(N.x + N.y)
chrom_stat <- rename(chrom_stat, Nm=N.x,Nf=N.y, Male.Hom=Obs.Hom.x, Fem.Hom=Obs.Hom.y)
chrom_stat$weight <- chrom_stat$Nf + chrom_stat$Nm

# Excluding unordered contigs
chrom_stat <- chrom_stat[grep("chr.*", chrom_stat$Chr),]

ggplot(chrom_stat, aes(x=BP, y=hom)) + facet_grid(~Chr, scales='free_x') + geom_smooth(method = 'lm')
