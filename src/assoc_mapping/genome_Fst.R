# This scripts is used to visualize the smoothed Fst value per loci along the genome.
# This is an exploratory analyse prior to association mapping.

library(dplyr)

in_path <- commandArgs(TrailingOnly=T)[1]
phi_stat <- read.csv('../../data/populations/d-25_r-75/batch_0.phistats_F-M.tsv',header=T,skip=2,sep='\t')

phi_stat <- tbl(phi_stat)
by_chrom <- phi_stat %>% group_by(Chr) %>% summarise(avg=mean(Fst.))
plot(by_chrom,type='p')

