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


hom_filt <- function(fam){
# Function for filtering loci homozygous in mother of given family
  mother <- as.character(indv$Name[indv$Family==fam & indv$Generation=='F3'])
  # Using mother's genotype from snps file
  snp_gz <- gzfile(paste0('../../data/sstacks/',fam,'/',mother,'.snps.tsv.gz'),"rt")
  snp_mother <- read.table(snp_gz); close(snp_gz)
  # keeping all loci where all SNPs are hom
  hom_loc <- snp_mother %>% group_by(V3) %>% filter(all(V5=='O'))
  hom_ID_local <- unique(hom_loc$V3)
  match_gz <- gzfile(paste0('../../data/sstacks/',fam,'/',mother,'.matches.tsv.gz'),"rt")
  match_mother <- read.table(match_gz); close(match_gz)
  hom_ID_global <- unique(match_mother$V3[match_mother$V4 %in% hom_ID_local])
  return(hom_ID_global)
}

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
# Customized operator for convenience

# Genome statistics
stat_path <- '../../data/populations/d-20_r-80/'
indv <- read.table('../../data/individuals', header=T)
# phi_path <- commandArgs(TrailingOnly=T)[1]
sum_stat <- data.frame()
for(fam in list.dirs(stat_path)[2:length(list.dirs(stat_path))]){  # Excluding first dir (parent)
  tmp_stat <- read.csv(paste0(fam,'/batch_0.sumstats.tsv'),header=T,skip=2,sep='\t')
  tmp_stat$fam <- rep(basename(fam))
  try(mother_hom <- hom_filt(basename(fam)))
  tmp_stat <- tmp_stat[tmp_stat$Locus.ID %not in% mother_hom,]
  sum_stat <- rbind(sum_stat, tmp_stat)
  mother_hom <- ""
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
