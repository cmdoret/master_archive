# This script parses the populations output vcf file to analyse and visualise coverage across the genome.
# Cyril Matthey-Doret
# 27.07.2017


library(ggplot2)
library(dplyr)
library(tidyr)

depth <- read.table("../../data/coverage/site_depth.txt", header=T)  # Loading data

# Excluding unordered contigs
depth <- depth[grep("chr*",depth$CHROM),]
depth$CHROM <- droplevels(depth$CHROM)

# Grouping SNPs every 1kb
depth$POS <- round(depth$POS, digits=-3)
depth <- depth %>%
  group_by(CHROM, POS, Family) %>%
  mutate(avg=mean(MEAN_DEPTH)) %>%
  select(-MEAN_DEPTH) %>%
  rename(MEAN_DEPTH=avg)

depth <- subset(depth, !duplicated(subset(
          depth, select=c(CHROM, POS, Family))))
  
# Spread families into columns
sp_depth <- spread(depth,Family,MEAN_DEPTH)
sp_depth <- sp_depth[rowSums(!is.na(sp_depth[,3:14]))>0,]

# Plotting coverage by family
ggplot(data=depth,aes(x=POS,y=MEAN_DEPTH,col=Family)) + 
  facet_grid(CHROM~., scales = 'free_x') + geom_point(size=0.3) + 
  geom_line(alpha=0.4) + theme_bw() + ggtitle("Coverage across genome, 1kb windows") + 
  xlab("Genomic position of SNP") + ylab("Mean coverage")


if(FALSE){
mergecov <- depth
ggplot(data=depth,(aes(x=POS,y=MEAN_DEPTH,col=Family))) + 
  facet_grid(~CHROM,scales = 'free_x') + geom_point(alpha=0.3)

slide <- data.frame()
for(chr in levels(mergecov$CHROM)){
  for(bp in unique(mergecov$POS[mergecov$CHROM==chr])){
    curr <- data.frame(chr,bp,round(mean(mergecov$MEAN_DEPTH[mergecov$POS==bp & mergecov$CHROM==chr]),3))
    slide <- rbind(slide,curr)
  }
}
colnames(slide) <- c("CHROM","POS","DEPTH")

ggplot(data=slide,aes(x=POS,y=DEPTH)) + geom_line() + facet_wrap(~CHROM,scales = 'free_x')
}
