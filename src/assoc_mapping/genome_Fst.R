# This scripts is used to visualize the smoothed Fst value per loci along the genome.
# This is an exploratory analyse prior to association mapping.

# Loading and formatting data

library(dplyr)
library(readr)
library(ggplot2)


# Chromosome sizes
scaffolds <- read.table("../../data/ref_genome/ordered_genome/merged.fasta.ann")
chrom_sizes <- data.frame(); chrom_names=c()
for(row in seq(1,nrow(scaffolds))){
  if(length(grep('chr',scaffolds[row,2]))){
    chrom_names <- append(chrom_names,as.character(scaffolds[row,2]))
    chrom_sizes <- rbind(chrom_sizes, scaffolds[row+1,c(1,2)])
  }
}
chrom_sizes <- cbind(chrom_names, chrom_sizes)
colnames(chrom_sizes) <- c("chrom","start","length")

# Genome statistics
phi_path <- '../../data/populations/d-20_r-80/'
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

highFST <- chrom %>%
  group_by(fam) %>%
  filter(Fst.==max(Smoothed.Fst.))


ggplot(data=chrom, aes(x=tot_BP, y=Smoothed.Fst.))+ geom_vline(data=chrom_sizes, aes(xintercept=start), col="blue", lty=2) + 
  geom_line() + facet_wrap(~fam,drop = F) + geom_text(data=highFST,aes(label=pos, y=Fst.+0.2), size=1.9)
