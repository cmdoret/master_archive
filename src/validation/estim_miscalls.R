# In this script, I verify the validity of CSD association hits using the number of heterozygous females at each candidate.
# Cyril Matthey-Doret
# 04.12.2017


#=== LOADING ===#
# Loading packages
pack <- c("ggplot2","dplyr", "readr", "viridis", "reshape2", "tibble")
lapply(pack, require, character.only = TRUE)
# Loading E/O/M genotype matrix
hom_path <- "data/assoc_mapping/grouped_keep_geno_EOM.tsv"
sum_stat <- read_tsv(file=hom_path, col_names=T)
# Loading list of individuals with family, generation and sex informations
indv <- read.table("data/individuals.tsv", header=T, sep='\t', stringsAsFactors = F)
# Loading list of SNPs with significant p-values for CSD association mapping
hits <- read.table("data/assoc_mapping/case_control/case_control_hits.tsv", header=T, stringsAsFactors = F)

#=== PROCESSING ===#
# Filtering most significant hit on each chromosome
best <-  hits %>%
  group_by(Chr) %>%
  filter( fisher > 3 & fisher == max(fisher))

# Subset best CSD SNPs
csd_stat <- sum_stat %>%
  right_join(., best, by=c("Chr","BP")) %>%
  select(2,3,one_of(indv$Name[indv$Sex=='F'])) %>%
  mutate_at(funs(ifelse(.=='E',1,0)),.vars = vars(-Chr,-BP))

#csd_stat <- aggregate(csd_stat, by=list(csd_stat$Chr),FUN=max)


# Convert each column to a byte containing N bits where N is the number of CSD loci
# E = het = 1; O = hom = 0
# Single loci encoding will be: 1, 2, 4, 8, ...
binary_vals <- 2^(0:(nrow(csd_stat)-1))
csd_stat[,3:dim(csd_stat)[2]] <- csd_stat[,3:dim(csd_stat)[2]] * binary_vals
csd_sum <- colSums(csd_stat[,3:dim(csd_stat)[2]])

id <- paste(csd_stat$Chr,csd_stat$BP,sep=':')
csd_mat <- matrix(nrow=nrow(csd_stat),ncol = nrow(csd_stat),dimnames = list(id,id))
# Looping over all pairwise combinations of loci
for ( i in 1:nrow(csd_stat) ) {
  for (j in 1:nrow(csd_stat) ) {
    if ( i != j){
      # On the diagonal: Heterozygous only at locus i (==j)
      csd_mat[i,j] <- length(csd_sum[csd_sum == 2^(i-1) + 2^(j-1) ])
    }else{
      # Position i,j in the matrix: heterozygous at loci i and j only
      csd_mat[i,j] <- length(csd_sum[csd_sum == 2^(i-1) ])
    }
  }
}
csd_plot <- melt(csd_mat)

ggplot(data=csd_plot,aes(Var1,Var2,fill=value)) + geom_tile() + 
  scale_fill_viridis() + labs(fill = "Females") + xlab("") + ylab("")

  
# Estimating mis-calling rate using portion of unlikely homo -> het transitions
# 2 METRICS:
#=== ASSUMING CORRECT CALLS IN MOTHERS ===# --> wrong_F4
# For all mother-homozygous SNPs, count the proportion of offspring calls that are heterozygous.
#=== "NOT ASSUMPTION", THRESHOLD INSTEAD ===# --> wrong_any
# Each SNP can have at most 1 wrong call. Either it is wrong in an offspring, or in the mother (where it was
# in fact homozygous). Small assumption: Unlikely to have genotype miscalls twice at the same SNP.

# Excluding chromosome and bp infos
snps <-sum_stat[,4:ncol(sum_stat)] %>% 
  # Exclude uninformative fixed homozygous sites (will speed up next step)
  filter(rowSums(.=='E')>0) %>%
  # filter only those where at least one mother is heterozygous
  filter_at(vars(one_of(indv$Name[indv$Generation=='F3'])), any_vars(. == 'O')) %>%
  # Transposing data -> snps are columns, samples are rows
  t(.) %>%
  data.frame(.,stringsAsFactors = F) %>%
  # Names used to merge with individuals to add Sex, generation and family infos
  rownames_to_column("Name") %>% 
  merge(indv,by="Name")


# Total number of offspring samples per family  
stat_snp <- snps %>%
  group_by(Family) %>% 
  summarise(samples=n()) %>%
  mutate(wrong_any = 0, wrong_F4 = 0, correct = 0)
for ( fam in unique(snps$Family)){
  
  tidy_snp <- snps %>%
    filter(Family==fam) %>%
    select(-Name, -Family, -Sex) %>%
    reshape2::melt(id.vars = "Generation")
  
  F3O <-  tidy_snp %>% 
    filter(value == "O" & Generation == "F3") %>% 
    select(variable)
  
  F4E <- tidy_snp %>%
    filter(value == "E" & Generation == "F4") %>% 
    select(variable)
  
  F4O <- tidy_snp %>%
    filter(value == "O" & Generation == "F4") %>% 
    select(variable)
  
  wrong_snps <- intersect(unique(F4E$variable),unique(F3O$variable))
  correct_snps <- intersect(unique(F4O$variable),unique(F3O$variable))
  
  wrong <- snps %>% 
    filter(Generation == "F4" & Family == fam) %>%
    select(wrong_snps)
  
  correct <- snps %>% 
    filter(Generation == "F4" & Family == fam) %>%
    select(correct_snps)
  
  stat_snp$wrong_F4[stat_snp$Family==fam] <- sum(apply(wrong, MARGIN = 2,function(x) sum(x == "E")))
  stat_snp$correct[stat_snp$Family==fam] <- sum(apply(correct, MARGIN = 2,function(x) sum(x == "O")))
  stat_snp$wrong_any[stat_snp$Family==fam] <- sum(apply(wrong, MARGIN = 2,function(x) any(x == "E")))

}

stat_snp <- filter(stat_snp, wrong_any > 0 | correct > 0)

plt_snp <- stat_snp %>% 
  mutate(wrong_F4 = wrong_F4/(wrong_F4+correct), 
         wrong_any = wrong_any/(wrong_any+correct)
         ) %>% 
  select(-correct, -samples) %>%
  melt(id.vars = "Family",
     variable.name = "Measurement",
     value.name = "Miscalls")

ggplot(data=plt_snp, aes(x=Family,fill=Measurement, y=Miscalls)) + 
  geom_bar(stat="identity", position="dodge") + theme_bw() + ggtitle("Proportion of wrong genotype calls") + 
  scale_x_discrete(labels=sprintf("%s: %i", Family, stat_snp$samples)) + 
  scale_fill_brewer(palette = "Paired", labels=c("P(F4:E > 0)","P(F4:E | F3:O)")) 
