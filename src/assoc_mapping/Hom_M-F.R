# This scripts is used to visualize the CSD-ness of SNPs along the genome.
# This is measured as [Obs.Hom(M)+Obs.Het(F)]/2
# 09.07.2017
# Cyril Matthey-Doret

# Loading and formatting data

library(dplyr)
library(readr)
library(ggplot2)

hom_path <- "../../data/assoc_mapping/prop_hom_fixed_sites.tsv"
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


sum_stat <- read_tsv(hom_path, col_names=T)
sum_stat <- sum_stat[sum_stat$N.Samples>0,]

# Computing CSD-ness
CSD_like <- sum_stat %>% mutate(CSD=((1-Prop.Hom.F)+Prop.Hom.M)/2)

# Keeping only chromosomes (removing contigs)
CSD_like <- CSD_like %>% arrange(Chr, BP)
chrom <- CSD_like[grep("chr.*", CSD_like$Chr),]
chrom$pos <- paste(chrom$Chr,chrom$BP, sep="_")
chrom$Chr <- as.factor(chrom$Chr)
chrom$Chr <- droplevels(chrom$Chr)

# Finding genomic position of SNPs
genomic_pos <- function(snp){
  # Computing total basepair position in genome as start(chrom)+bp
  gen_pos <- as.numeric(snp["BP"]) + chrom_sizes[chrom_sizes$chr==snp["Chr"],"start"]
  return(gen_pos)
}
chrom$tot_BP <-apply(X = chrom,MARGIN = 1, FUN=genomic_pos)

# Proportion of offspring homozygous at loci heterozygous in mother
ggplot(data=chrom,aes(x=BP,y=CSD)) + geom_point() +facet_grid(~Chr, scales = "free_x") 

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
  facet_wrap(~Family,drop = F) + #geom_text(data=highCSD,aes(label=pos, y=CSD), size=1.9)+
  geom_vline(data=chrom_sizes, aes(xintercept=start), col="blue", lty=2) +
  ylab("Prop.CSD")

ggplot(data=chrom, aes(x=tot_BP, y=CSD, col=Family))+ geom_line() + 
  #geom_text(data=highCSD,aes(label=pos, y=CSD), size=1.9)+
  geom_vline(data=chrom_sizes, aes(xintercept=start), col="blue", lty=2) +
  ylab("Prop.CSD")

uniq_CSD = top_CSD[order(top_CSD[,'tot_BP'],-top_CSD[,'CSD']),]
uniq_CSD$tot_BP <- round(uniq_CSD$tot_BP,digits=-3)
uniq_CSD$BP <- round(uniq_CSD$BP,digits=-3)
uniq_CSD = uniq_CSD[!duplicated(uniq_CSD$tot_BP),]

out_CSD <- uniq_CSD[,c("Locus.ID","Chr","BP","tot_BP","Nf","Nm","CSD")]
write.csv(out_CSD,"../../data/assoc_mapping/CSD_hits.csv",row.names=F,quote=F)
