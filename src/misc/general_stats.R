# General statistics, per individual and per locus to have an overview of the data.
# Cyril Matthey-Doret
# 14.07.2017

library(ggplot2); library(ggjoy)

# Per individual

indv <- read.table("../../data/ploidy/thresholds/fixed",header = T)
indv$lib <- "Mothers"
for(i in c("6","7","7b","10","10b")){
  tmp_lib <- read.table(paste0("../../data/barcodes/barcodes_radwasp", i))
  indv$lib[indv$Name %in% tmp_lib[,2]] <- i
}

ggplot(data=indv, aes(x=N_SITES, y=lib, height=..density..)) + 
  geom_joy(stat = "binline", bins=20, scale=0.95, draw_baseline=T)+
  ylab("Library")

