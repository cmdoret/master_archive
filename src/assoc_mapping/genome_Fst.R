# This scripts is used to visualize the smoothed Fst value per loci along the genome.
# This is an exploratory analyse prior to association mapping.

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
phi_path <- '../../data/populations/d-20_r-80/'
# phi_path <- '../../data/populations/haplo_d-20_r-80/'
# phi_path <- commandArgs(TrailingOnly=T)[1]
phi_stat <- data.frame()
for(fam in list.dirs(phi_path)[2:length(list.dirs(phi_path))]){  # Excluding first dir (parent)
tmp_phi <- read.csv(paste0(fam,'/batch_0.phistats_F-M.tsv'),header=T,skip=2,sep='\t')
tmp_phi$fam <- rep(basename(fam))
phi_stat <- rbind(phi_stat, tmp_phi)
}

# Keeping only chromosomes (removing contigs)
phi_stat <- phi_stat %>% arrange(Chr, BP)
chrom <- phi_stat[grep("chr.*", phi_stat$Chr),]
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
  summarise(avg=mean(Smoothed.Fst.), BP=mean(tot_BP)) %>%
  arrange(BP)

plot(compact_chrom$BP, compact_chrom$avg, type="l", ylab='Fst', xlab="genomic position")
abline(v=chrom_sizes$start, lty=2,col="blue")
text(x = chrom_sizes$mid,y=-0.2,labels = chrom_sizes$chrom)
hit_ID <-compact_chrom$Locus.ID[order(compact_chrom$avg,decreasing = T)][1:3]
text(x = chrom$tot_BP[chrom$Locus.ID %in% hit_ID],
     y= compact_chrom$avg[compact_chrom$Locus.ID %in% hit_ID],
     labels = chrom$pos[chrom$Locus.ID %in% hit_ID])

highFST <- chrom %>%
  group_by(fam) %>%
  filter(Smoothed.Fst.==max(Smoothed.Fst.))

ggplot(data=chrom, aes(x=tot_BP, y=Smoothed.Fst.))+ geom_vline(data=chrom_sizes, aes(xintercept=start), col="blue", lty=2) + 
  geom_line() + facet_wrap(~fam,drop = F) + geom_text(data=highFST,aes(label=pos, y=Fst.+0.2), size=1.9)

top_Fst<- chrom[chrom$Fst.>=0.7,]
hist(top_Fst$tot_BP,breaks=100, main="Top CSD candidates", xlab="Genomic position", ylab="N hits >= 0.8", col="grey")
abline(v=chrom_sizes$start, lty=2,col="blue")

uniq_Fst = top_Fst[order(top_Fst[,'tot_BP'],-top_Fst[,'Fst.']),]
uniq_Fst$tot_BP <- round(uniq_Fst$tot_BP,digits=-3)
uniq_Fst$BP <- round(uniq_Fst$BP,digits=-3)
uniq_Fst = uniq_Fst[!duplicated(uniq_Fst$tot_BP),]

out_Fst <- uniq_Fst[,c("Locus.ID","Chr","BP","tot_BP")]
write.csv(out_Fst,"../../data/assoc_mapping/Fst_hits.csv",row.names=F,quote=F)
