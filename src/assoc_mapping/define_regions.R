# Script used to define the size of CSD candidate regions using stats of RADseq SNPs.
# Cyril Matthey-Doret
# 08.03.2018

# LOAD
# Packages
packs <- c("ggplot2","dplyr", "gridExtra", "readr")
lapply(packs, require, character.only=T)
# Data
snps <- read.table("data/assoc_mapping/case_control/case_control_all.tsv", header=T)
# User defined p-value significance threshold, in -log10(pvalue)
thresh <- 3



# PROCESS
# Mark significant SNPs as part of CSD region
snps <- snps %>%
  mutate(csd_locus = ifelse(fisher>thresh, yes=1, no=0))

region <- 1
# Iterate over all snps
for(row in 1:nrow(snps)){
  # If snps is significant
  if(snps$csd_locus[row]){
    # Assign it a region number
  if(snps$csd_locus[row+1]){
    snps$csd_locus[row] <- region
  }
  # When region ends set new region number
  else{region <- region + 1}
  }
}

region <- 1
for(row in 1:nrow(snps)){
  # If SNP is significant
  if(snps$csd_locus[row]){
    # And next one is also significant, consider region
    snps$csd_locus[row] <- ifelse(snps$csd_locus[row+1], yes=region, no=0)
    # end of region, 
    if(snps$csd_locus[row-1] & !snps$csd_locus[row+1]) {
      # change region number for next one
      region <- region + 1
    }
  }
}

# VISUALISE
ggplot(data=snps, aes(x=BP, y=fisher)) + 
  geom_point() +
  geom_point(data=snps[snps$csd_locus>thresh,], 
             aes(col=as.factor(csd_locus)), 
             size=3) +
  geom_line(data=snps[snps$csd_locus>thresh,], 
             aes(col=as.factor(csd_locus))) +
  facet_grid(~Chr)
