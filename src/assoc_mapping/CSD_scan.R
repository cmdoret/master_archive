# This scripts is used to scan SNPs along the genome for CSD-ness.
# This is measured as [Obs.Hom(M)+Obs.Het(F)]/2
# 09.07.2017
# Cyril Matthey-Doret

# Loading and formatting data

library(dplyr)
library(ggplot2)
library(readr)
#hom_path <- "../../data/assoc_mapping/grouped_prop_hom_fixed_sites.tsv"
hom_path <- commandArgs(trailingOnly = T)[1]  # Path to process_genomic.py output
ref_genome <- commandArgs(trailingOnly = T)[2]  # Path to reference genome ann file
out_path <- commandArgs(trailingOnly = T)[3] # Output path
thresh <- commandArgs(trailingOnly = T)[4]  # Hard filter for CSD hits
scaffolds <- read.table(ref_genome, stringsAsFactors = F)
indv <- read.table('data/individuals.tsv', header=T)

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
#sum_stat <- read.table(hom_path, header=T, na.strings='NA', sep='\t')
sum_stat <- read_tsv(file=hom_path, col_names=T, col_types = "iciddddddc")
sum_stat <- sum_stat[sum_stat$N.Samples>0,]
sum_stat <- sum_stat[!is.na(sum_stat$Prop.Hom),]
# Computing CSD-ness
CSD_like <- sum_stat %>% mutate(CSD=((1-Prop.Hom.F)+Prop.Hom.M)/2)

# Keeping only chromosomes (removing unordered contigs)
CSD_like <- CSD_like %>% arrange(Chr, BP)
chrom <- CSD_like[grep("chr.*", CSD_like$Chr),]
chrom$pos <- paste(chrom$Chr,chrom$BP, sep="_")
chrom$Chr <- as.factor(chrom$Chr)
chrom$Chr <- droplevels(chrom$Chr)
chrom <- chrom[chrom$N.Females>0 & chrom$N.Males>0,]

# Finding genomic position of SNPs
genomic_pos <- function(snp){
  # Computing total basepair position in genome as start(chrom)+bp
  gen_pos <- as.numeric(snp["BP"]) + chrom_sizes[chrom_sizes$chr==snp["Chr"],"start"]
  return(gen_pos)
}
chrom$tot_BP <-apply(X = chrom,MARGIN = 1, FUN=genomic_pos)

# Overview: All families on same plot with smoothing
# ggplot(data=chrom[chrom$N.Females>2 & chrom$N.Males>2,],aes(x=BP,y=CSD, col=Family, weight=N.Samples)) +facet_grid(~Chr, scales = "free_x") + 
#   stat_smooth(span=0.75)
# # Pooling all families together (averaging each SNP)
# # not getting anything from it, different het CSD loci in different families.
# compact_chrom <- chrom[chrom$Family!='J',] %>%
#   group_by(Locus.ID) %>%
#   summarise(avg=mean(CSD), indv=length(unique(Family)), BP=mean(tot_BP)) %>%
#   arrange(BP)
# plot(compact_chrom$BP, compact_chrom$avg, type="l", ylab='prop. CSD', xlab="genomic position",
#      main="Proportion of CSD individuals, averaged across families")
# mod <- loess(data=compact_chrom, formula=avg~BP, weights = indv, span = 0.15, degree = 2, model=T)
# points(mod$x[order(mod$x)],mod$fitted[order(mod$x)], type='l',lwd=1.3, col='red')
# abline(v=chrom_sizes$start, lty=2,col="blue")
# text(x = chrom_sizes$mid,y=0.2,labels = chrom_sizes$chrom)

pdf(paste0(out_path, "/plots/scan_CSD_hits_", thresh, ".pdf"))
top_CSD <- chrom[chrom$CSD>=thresh,]
hist(top_CSD$tot_BP,breaks=100, main="Top CSD candidates", xlab="Genomic position", ylab=paste0("N hits >= ", thresh), col="grey")
abline(v=chrom_sizes$start, lty=2,col="blue")
text(x = chrom_sizes$mid,y=17,labels = chrom_sizes$chrom)

pdf(paste0(out_path, "/plots/scan_CSD_per_fam.pdf"))
ggplot(data=chrom, aes(x=tot_BP, y=CSD))+ geom_line() +
  facet_wrap(~Family,drop = F) + #geom_text(data=highCSD,aes(label=pos, y=CSD), size=1.9)+
  geom_vline(data=chrom_sizes, aes(xintercept=start), col="blue", lty=2) +
  ylab("Prop.CSD")
dev.off()
# ggplot(data=chrom, aes(x=tot_BP, y=CSD, col=Family))+ geom_line() + 
#   #geom_text(data=highCSD,aes(label=pos, y=CSD), size=1.9)+
#   geom_vline(data=chrom_sizes, aes(xintercept=start), col="blue", lty=2) +
#   ylab("Prop.CSD")

#uniq_CSD = top_CSD[order(top_CSD[,'tot_BP'],-top_CSD[,'CSD']),]
#uniq_CSD$tot_BP <- round(uniq_CSD$tot_BP,digits=-3)
#uniq_CSD$BP <- round(uniq_CSD$BP,digits=-3)
#uniq_CSD = uniq_CSD[!duplicated(uniq_CSD$tot_BP),]

out_CSD <- top_CSD
write.table(out_CSD,paste0(out_path, "/hits/CSD_scan_hits_", thresh, ".tsv"),row.names=F,quote=F)
