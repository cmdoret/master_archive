# This scripts is used to visualize the smoothed Fst value per loci along the genome.
# This is an exploratory analyse prior to association mapping.

library(dplyr)
library(readr)

in_path <- commandArgs(TrailingOnly=T)[1]
in_path <- '../../data/populations/d-20_r-75/'

par(mfrow=c(2,6))
phi_stat <- data.frame()

for(fam in list.dirs(in_path)[2:length(list.dirs(in_path))]){  # Excluding first dir (parent)
tmp_phi <- read.csv(paste0(fam,'/batch_0.phistats_F-M.tsv'),header=T,skip=2,sep='\t')
tmp_phi$fam <- rep(basename(fam))
phi_stat <- rbind(phi_stat, tmp_phi)
}

phi_stat <- phi_stat %>% arrange(Chr, BP)
by_chrom <- phi_stat %>% group_by(Chr) %>% summarise(avg=mean(Fst.))
#plot(by_chrom,type='p')

chrom <- phi_stat[grep("chr.*", phi_stat$Chr),]
chrom$pos = paste(chrom$Chr,chrom$BP, sep="_")

plot(rownames(chrom), chrom$Fst., type="l", ylab='Fst')
highFST <- chrom[chrom$Fst.==max(chrom$Fst.),]
text(x = as.numeric(rownames(highFST)),y = highFST$Fst.,labels = highFST$pos)
chrom$Chr=droplevels(chrom$Chr)
for(i in levels(chrom$Chr)){
  tmp_chrom <- chrom[chrom$Chr==i,]
  abline(v=min(as.numeric(rownames(tmp_chrom))), col='blue')
}
