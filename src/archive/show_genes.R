# Just a script to produce genome browser-like visualization from annotations

library(ggplot2)
library(dplyr)

annot <- read.table('data/annotations/csd_gwas_hit_annot.tsv', sep='\t', header=T)

hits_df <- read.table('data/assoc_mapping/case_control/case_control_hits.tsv', sep='\t', header=T)
hits_df <- hits_df %>% rename(chrom=Chr)
annot <- annot %>% 
  group_by(chrom) %>%
  mutate(rstart=start)
ggplot(annot) + geom_rect(aes(xmin = start, xmax = end, ymin = -0.3, ymax = -0.1)) + facet_grid(~chrom, scales='free_x') + 
  geom_hline(aes(yintercept=-0.2)) + geom_point(data=hits_df,aes(x=BP,y=-0.2),col='red') + theme_bw() + 
  ylab("Annotations per 10kb") + geom_histogram(data=annot,aes(x=start), binwidth=10000)
                                                                                                         