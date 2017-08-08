# This scripts is used to visualize the CSD-ness of SNPs along the genome.
# This is measured as [Obs.Hom(M)+Obs.Het(F)]/2
# 09.07.2017
# Cyril Matthey-Doret

# Loading and formatting data

library(dplyr)
library(readr)
library(ggplot2)

stat_path <- '../../data/populations/d-5_r-10/'
indv <- read.table('../../data/individuals', header=T)
scaffolds <- read.table("../../data/ref_genome/ordered_genome/merged.fasta.ann", stringsAsFactors = F)
grouped <- "T"
# Summarizing chromosome sizes

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
chrom_sizes$mid <- chrom_sizes$start + (chrom_sizes$length/2)  # Middle position in chrom. used for plotting


if(grouped == 'F'){
  sum_stat <- data.frame()
  for(fam in list.dirs(stat_path)[2:length(list.dirs(stat_path))]){  # Excluding first dir (parent)
    tmp_stat <- read.csv(paste0(fam,'/batch_0.sumstats.tsv'),header=T,skip=2,sep='\t')
    tmp_stat$fam <- rep(basename(fam))
    sum_stat <- rbind(sum_stat, tmp_stat)
    mother_hom <- NULL
  }
  # removing SNPs homozygous in mother
  rm_snps<- sum_stat %>% 
    group_by(fam,Locus.ID,Col) %>%  # Each group contains one SNP (both male and female pop)
    summarise(hom_mot = sum(Pop.ID=='F' & Q.Nuc=='-'))  # If all females are homozygous -> mother is homozygous
  
  rm_snps <- rm_snps[rm_snps$hom_mot==0, 1:3]
  sum_stat <- merge(sum_stat,rm_snps,by=c("Locus.ID","Col","fam"))
} else{
  sum_stat <- read.csv(paste0(stat_path, '/batch_0.sumstats.tsv'), header=T, skip=2, sep='\t')
  sum_stat$fam <- 'all'
}

# Computing CSD-ness
male_stat <- sum_stat[sum_stat$Pop.ID=="M",c('Chr','BP','Obs.Hom', 'fam', 'N',"Locus.ID")]
female_stat <- sum_stat[sum_stat$Pop.ID=="F",c('Chr','BP','Obs.Het', 'fam', 'N')]
CSD_like <- merge(male_stat,female_stat,by=c('Chr','BP','fam'))
#CSD_like$CSD <- (CSD_like$Obs.Hom+CSD_like$Obs.Het)*(CSD_like$N.x+CSD_like$N.y)/2
CSD_like$CSD <- (CSD_like$Obs.Hom+CSD_like$Obs.Het)/2
CSD_like <- rename(CSD_like, Nm=N.x,Nf=N.y, Male.Hom=Obs.Hom, Fem.Het=Obs.Het)
#CSD_like <- CSD_like[CSD_like$Nf>1 & CSD_like$Nm>1,]


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

# Proportion of offspring homozygous at loci heterozygous in mother
ggplot(data=chrom,aes(x=BP,y=((1-Fem.Het)+Male.Hom))) + geom_point() +facet_grid(~Chr, scales = "free_x") 

compact_chrom <- chrom %>% 
  group_by(Locus.ID) %>%
  summarise(avg=mean(CSD), BP=mean(tot_BP)) %>%
  arrange(BP)

plot(compact_chrom$BP, compact_chrom$avg, type="l", ylab='prop. CSD', xlab="genomic position", 
     main="Proportion of CSD individuals, averaged across families")
abline(v=chrom_sizes$start, lty=2,col="blue")
text(x = chrom_sizes$mid,y=0.2,labels = chrom_sizes$chrom)

highCSD <- chrom %>%
  group_by(fam) %>%
  filter(CSD==max(CSD))

top_CSD <- chrom[chrom$CSD>=0.8,]
hist(top_CSD$tot_BP,breaks=100, main="Top CSD candidates", xlab="Genomic position", ylab="N hits >= 0.8", col="grey")
abline(v=chrom_sizes$start, lty=2,col="blue")
text(x = chrom_sizes$mid,y=17,labels = chrom_sizes$chrom)

ggplot(data=chrom, aes(x=tot_BP, y=CSD))+ geom_line() +
  facet_wrap(~fam,drop = F) + #geom_text(data=highCSD,aes(label=pos, y=CSD), size=1.9)+
  geom_vline(data=chrom_sizes, aes(xintercept=start), col="blue", lty=2) +
  ylab("Prop.CSD")

ggplot(data=chrom, aes(x=tot_BP, y=CSD, col=fam))+ geom_line() + 
  #geom_text(data=highCSD,aes(label=pos, y=CSD), size=1.9)+
  geom_vline(data=chrom_sizes, aes(xintercept=start), col="blue", lty=2) +
  ylab("Prop.CSD")

uniq_CSD = top_CSD[order(top_CSD[,'tot_BP'],-top_CSD[,'CSD']),]
uniq_CSD$tot_BP <- round(uniq_CSD$tot_BP,digits=-3)
uniq_CSD$BP <- round(uniq_CSD$BP,digits=-3)
uniq_CSD = uniq_CSD[!duplicated(uniq_CSD$tot_BP),]

out_CSD <- uniq_CSD[,c("Locus.ID","Chr","BP","tot_BP","Nf","Nm","CSD")]
write.csv(out_CSD,"../../data/assoc_mapping/CSD_hits.csv",row.names=F,quote=F)
