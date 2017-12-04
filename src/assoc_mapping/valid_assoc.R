# In this script, I verify the validity of CSD association hits using the number of heterozygous females at each candidate.
# Cyril Matthey-Doret
# 04.12.2017


#=== LOADING ===#
# Loading packages
pack <- c("ggplot2","dplyr", "readr")
lapply(pack, require, character.only = TRUE)
# Loading genotype matrix
hom_path <- "data/assoc_mapping/grouped_keep_geno_EOM.tsv"
sum_stat <- read_tsv(file=hom_path, col_names=T)
# Loading list of individuals
indv <- read.table("data/individuals.tsv", header=T, sep='\t', stringsAsFactors = F)
# Loading list of significant p-values
hits <- read.table("data/assoc_mapping/case_control/case_control_hits.tsv", header=T, stringsAsFactors = F)

#=== PROCESSING ===#
# Filtering significant hits
best <-  hits %>%
  group_by(Chr) %>%
  filter(fisher == max(fisher) & fisher > 3)

# Subset CSD SNPs
csd_stat <- sum_stat %>%
  right_join(., best, by=c("Chr","BP")) %>%
  select(2,3,one_of(indv$Name[indv$Sex=='F'])) %>%
  mutate_at(funs(ifelse(.=='E',1,0)),.vars = vars(-Chr,-BP))

# Convert each column to a byte containing N bits where N is the number of CSD loci
# E = het = 1; O = hom = 0
# Single loci encoding will be: 1, 2, 4, 8, ...
binary_vals <- 2^(0:(nrow(csd_stat)-1))
csd_stat[,3:dim(csd_stat)[2]] <- csd_stat[,3:dim(csd_stat)[2]] * binary_vals
csd_sum <- colSums(csd_stat[,3:dim(csd_stat)[2]])
csd_num <- csd_stat[,1:2]
csd_mat <- matrix(nrow=nrow(csd_stat),ncol = nrow(csd_stat))
