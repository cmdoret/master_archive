# Analysis of RADseq SNPs using PI values from WGS data.
# Cyril Matthey-Doret
# 15.03.2018

library("ggplot2", "viridis")
wgs_pi <- read.table("data/wgs_wild_mothers/RAD_sites/rad2wgs_pi.vcf.sites.pi", header=T)
rad_snps <- read.table("data/assoc_mapping/case_control/case_control_all.tsv", header=T)
colnames(wgs_pi)[1:2] <- c("Chr","BP")

snps <- merge(x = wgs_pi, 
              y = rad_snps, 
              by = c("Chr", "BP"), 
              all = F)
plot(snps$fisher, snps$PI)

ggplot(snps, aes(x=BP, y=fisher, col=PI)) + 
  geom_point() +
  facet_grid(~Chr) +
  viridis::scale_color_viridis()
