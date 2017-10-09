# Simple R scripts counting the number of lines (SNPs) belonging to each
# family in a simple CSV list of SNPs.
# Cyril Matthey-Doret
# 03.07.2017

load_packages <- function(){all(require(readr), require(dplyr))}

if(!suppressMessages(require(tidyverse))){
  stopifnot(load_packages())
}

# SNP_path <- "../../data/SNP_lists/m2_hom_mother.txt"
if(length(commandArgs(trailingOnly=TRUE))==2){
  SNP_path <- commandArgs(trailingOnly = TRUE)[1]
  out_name <- commandArgs(trailingOnly = TRUE)[2]
  print(paste0("Processing SNPs list ", basename(SNP_path), " ..."))
} else{
  print("Wrong number of arguments, need inpput path and output filename. Exiting now!")
  quit()
}

SNPs <-  read_csv(SNP_path)

SNPs %>% 
  group_by(family) %>%
  summarise("N SNPs"=length(family)) %>%
  write_csv(path = paste0("reports/lab_book/cleaning_genomic_data/", out_name), col_names = T)
