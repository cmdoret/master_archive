# Just computing some summary statistics over the reference genome to include them in lab book.
# Cyril Matthey-Doret
# 11.06.2017

library(tidyverse)

# Path to reference genome
ref_path <- commandArgs(trailingOnly = T)[1]
summary_path <- "src/misc/ref_contigs.txt"

# Extracting contig headers
system(paste('cat', ref_path, '| grep ">" >', summary_path,sep=' '))
ref <- read_delim(file = summary_path, delim = ' ', col_names = F)[,1:2]
colnames(ref) <- c("contig","len")
ref$len %<>% 
  gsub(pattern = "len=([0-9]*)", replacement = "\\1") %>%
  as.numeric(.)

half_genome <- sum(ref$len)/2
ord_contigs <- ref %>% arrange(desc(len))
iter_N50 <- 0
for(n in 1:nrow(ord_contigs)){
  row <- ord_contigs[n,]
  iter_N50 <- iter_N50 + row$len
  if(iter_N50>=half_genome){
    n50=row$len
    break
  }
}

sum_stats <- data.frame(Statistics=c("Assembly length (Mbp)", "Largest contig (Mbp)", "Mean contig size (kbp)",
                                     "Median contig size (kbp)", "N50 (kbp)", "Number of contigs"),
                        Values=c(round(sum(ref$len)/1000000, 1), round(max(ref$len)/1000000,1), 
                                 round(mean(ref$len)/1000, 1), round(median(ref$len)/1000, 1), round(n50/1000, 1), length(ref$len)))

write.csv(sum_stats, file="reports/lab_book/ref_genome_stats.csv",row.names = F,quote = F)