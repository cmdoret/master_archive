# This script is used to visualize offspring proportions per female.
# It takes an input list of individuals with sex, family, generation 
# and ploidy informations.
# Cyril Matthey-Doret
# 19.06.2017

load_packages <- function(){all(require(dplyr), require(tidyr))}

if(!suppressMessages(require(tidyverse))){
  stopifnot(load_packages())
}

# in_path <- "../../data/ploidy/thresholds/m2"
in_path <- commandArgs(trailingOnly = T)[1]  # Taking ploidy list as CL argument
in_table <- read.table(in_path, header=T)



sum_table <- in_table %>%
  filter(Generation=='F4') %>%  # Excluding mothers (F3 generation)
  select(Sex, Family, Ploidy) %>%  # Subsetting only factors of interest
  group_by(Family) %>%  
  unite(identity, c(Family, Sex, Ploidy), remove=T, sep="_") %>%  # uniting 3 columns into 1 containing 3 informations
  table(.)  # Getting counts for each combination (all individuals together)

extr_ID <- function(family, sex, ploidy){
  # Function for extracting counts of haploid males/diploid males/females per family
  ID_count <- sum_table[names(sum_table)==paste(family, sex, ploidy, sep='_')]
  return(ifelse(test = length(ID_count)>0, yes = ID_count, no = 0))
}

store_fam_ploidy <- as.data.frame(matrix(nrow=length(levels(in_table$Family)), ncol = 4))
colnames(store_fam_ploidy) <- c('Family', 'diploid_son', 'haploid_son', 'daughter')
rownames(store_fam_ploidy) <- as.character(levels(in_table$Family))

pdf(paste0('reports/lab_book/ploidy_per_fam/', basename(in_path), ".pdf"))
par(mfrow=c(2, 6))
for(i in LETTERS){  # Producing a piechart for each family
  if(i %in% levels(in_table$Family)){
    fam_count <- rep(0, 3)
    fam_count[1] <- extr_ID(i, 'M', 'D')
    fam_count[2] <- extr_ID(i, 'M', 'H')
    fam_count[3] <- extr_ID(i, 'F', 'D')
    tmp_table <- data.frame(identity = c('Diploid Son', 'Haploid Son', 'Daughter'), 
                            count = fam_count)
    store_fam_ploidy[i,] <- append(i,fam_count)
    pie(x = tmp_table$count, labels = tmp_table$identity, main=i)
  }
}
dev.off()
write.table(store_fam_ploidy, file = paste0('reports/lab_book/ploidy_per_fam/', basename(in_path)), row.names = F)
