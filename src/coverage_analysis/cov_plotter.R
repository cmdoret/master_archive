# This script parses the populations output vcf file to analyse and visualise coverage across the genome.
# Cyril Matthey-Doret
# 27.07.2017


library(ggplot2)
library(dplyr)

depth <- read.table("../../data/coverage/site_depth.txt", header=T)
depth <- depth[grep("chr*",depth$CHROM),]
depth$CHROM <- droplevels(depth$CHROM)

mergecov <- depth
mergecov$POS <- round(mergecov$POS,digits=-4)

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

