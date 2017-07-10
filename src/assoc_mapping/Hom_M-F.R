# This scripts is used to visualize the CSD-ness of SNPs along the genome.
# This is measured as [Obs.Hom(M)+Obs.Het(F)]/2
# 09.07.2017
# Cyril Matthey-Doret

# Loading and formatting data

library(dplyr)
library(readr)
library(ggplot2)


# Chromosome sizes
scaffolds <- read.table("../../data/ref_genome/ordered_genome/merged.fasta.ann", stringsAsFactors = F)
chrom_sizes <- data.frame(); chrom_names=c()
for(row in seq(1,nrow(scaffolds))){
  if(length(grep('chr',scaffolds[row,2]))){
    chrom_names <- append(chrom_names,as.character(scaffolds[row,2]))
    chrom_sizes <- rbind(chrom_sizes, scaffolds[row+1,c(1,2)])
  }
}
chrom_sizes <- cbind(chrom_names, chrom_sizes)
rownames(chrom_sizes) <- NULL
colnames(chrom_sizes) <- c("chrom","start","length")
chrom_sizes$length <- as.numeric(chrom_sizes$length)
chrom_sizes$mid <- chrom_sizes$start + (chrom_sizes$length/2)

# Genome statistics
stat_path <- '../../data/populations/d-20_r-80/'
# phi_path <- commandArgs(TrailingOnly=T)[1]
sum_stat <- data.frame()
for(fam in list.dirs(stat_path)[2:length(list.dirs(stat_path))]){  # Excluding first dir (parent)
  tmp_stat <- read.csv(paste0(fam,'/batch_0.sumstats.tsv'),header=T,skip=2,sep='\t')
  tmp_stat$fam <- rep(basename(fam))
  sum_stat <- rbind(sum_stat, tmp_stat)
}


# Computing CSD-ness
male_stat <- sum_stat[sum_stat$Pop.ID=="M",c('Chr','BP','Obs.Hom', 'fam', 'N',"Locus.ID")]
female_stat <- sum_stat[sum_stat$Pop.ID=="F",c('Chr','BP','Obs.Het', 'fam', 'N')]
CSD_like <- merge(male_stat,female_stat,by=c('Chr','BP','fam'))
#CSD_like$CSD <- (CSD_like$Obs.Hom+CSD_like$Obs.Het)*(CSD_like$N.x+CSD_like$N.y)/2
CSD_like$CSD <- (CSD_like$Obs.Hom+CSD_like$Obs.Het)/2
CSD_like <- rename(CSD_like, Nm=N.x,Nf=N.y, Male.Hom=Obs.Hom, Fem.Het=Obs.Het)
CSD_like <- CSD_like[CSD_like$Nf>1 & CSD_like$Nm>1,]


# Keeping only chromosomes (removing contigs)
CSD_like <- CSD_like %>% arrange(Chr, BP)
chrom <- CSD_like[grep("chr.*", CSD_like$Chr),]
chrom$pos = paste(chrom$Chr,chrom$BP, sep="_")
chrom$Chr=droplevels(chrom$Chr)

# Finding genomic position of SNPs
genomic_pos <- function(snp){
  # Computing total basepair position in genome as start(chrom)+bp
  gen_pos <- as.numeric(snp["BP"]) + chrom_sizes[chrom_sizes$chr==snp["Chr"],"start"]
  return(gen_pos)
}
chrom$tot_BP <-apply(X = chrom,MARGIN = 1, FUN=genomic_pos)

compact_chrom <- chrom %>% 
  group_by(Locus.ID) %>%
  summarise(avg=mean(CSD), BP=mean(tot_BP)) %>%
  arrange(BP)


plot(compact_chrom$BP, compact_chrom$avg, type="l", ylab='prop. CSD', xlab="genomic position", 
     main="Proportion of CSD individuals, averaged across families")
abline(v=chrom_sizes$start, lty=2,col="blue")
text(x = chrom_sizes$mid,y=0,labels = chrom_sizes$chrom)

highCSD <- chrom %>%
  group_by(fam) %>%
  filter(CSD==max(CSD))

top_CSD <- chrom[chrom$CSD>=0.8,]
hist(top_CSD$tot_BP,breaks=100, main="Top CSD candidates", xlab="Genomic position", ylab="N hits >= 0.8", col="grey")
abline(v=chrom_sizes$start, lty=2,col="blue")
text(x = chrom_sizes$mid,y=17,labels = chrom_sizes$chrom)

ggplot(data=chrom, aes(x=tot_BP, y=CSD))+ geom_line() +
  facet_wrap(~fam,drop = F) + #geom_text(data=highCSD,aes(label=pos, y=CSD), size=1.9)+
  geom_vline(data=chrom_sizes, aes(xintercept=start), col="blue", lty=2) 

ggplot(data=chrom, aes(x=tot_BP, y=CSD, col=fam))+ geom_line() + 
  #geom_text(data=highCSD,aes(label=pos, y=CSD), size=1.9)+
  geom_vline(data=chrom_sizes, aes(xintercept=start), col="blue", lty=2) 

