# Just a script to produce genome browser-like visualization from annotations

library(ggplot2)
library(dplyr)

annot <- read.table('data/annotations/csd_gwas_hit_annot.tsv', sep='\t', header=T)

hits_df <- data.frame(chrom=c("tig00000010","tig00000010","tig00000010",
                              "tig00001764","tig00000797","tig00002024"),
                      BP=c(500869,495741,495732,10500,17843,132706))
annot <- annot %>% 
  group_by(chrom) %>%
  mutate(rstart=start
ggplot(annot) + geom_rect(aes(xmin = start, xmax = end, ymin = -1, ymax = 1)) + facet_grid(~chrom, scales='free_x') + 
  ylim(c(-10,10)) + geom_hline(aes(yintercept=0)) + geom_point(data=hits_df,aes(x=BP,y=0),col='red')
