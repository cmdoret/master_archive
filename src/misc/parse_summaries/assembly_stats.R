# Just computing some summary statistics over the reference genome to include them in lab book.
# Cyril Matthey-Doret
# 11.06.2017

load_packages <- function(){all(require(readr), require(dplyr))}

if(!suppressMessages(require(tidyverse))){
  stopifnot(load_packages())
}

# Path to reference genome
# ref_path <- "../../data/ref_genome/ordered_genome/merged.fasta"
# summary_path <- "../../src/misc/ref_ordered_contigs.txt"

ref_path <- commandArgs(trailingOnly = T)[1]
summary_path <- "src/misc/ref_ordered_contigs.txt"

# Extracting contig headers and lengths using awk one-liner:
# Note: $0 ~ ">" : ">" appears in the record
# cat ref_path | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > summary_path
system(paste("cat", ref_path, 
             "| awk '$0 ~ \">\" {print c; c=0;printf substr($0,2,100) \"\t\"; } $0 !~ \">\" {c+=length($0);} END { print c; }' > ", 
             summary_path, sep=" "))

ref <- read_delim(file = summary_path, delim = '\t', col_names = F)[,1:2]
colnames(ref) <- c("contig","len")

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



sum_stats <- data.frame(Statistics=c("Assembly length (Mbp)", "Largest scaffold (Mbp)", "Mean scaffold size (kbp)",
                                     "Median scaffold size (kbp)", "N50 (kbp)", "Number of scaffolds"),
                        Values=c(round(sum(ref$len)/1000000, 1), round(max(ref$len)/1000000,1), 
                                 round(mean(ref$len)/1000, 1), round(median(ref$len)/1000, 1), round(n50/1000, 1), length(ref$len)))

write.csv(sum_stats, file="reports/lab_book/ref_genome_stats.csv",row.names = F,quote = F)
