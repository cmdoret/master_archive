# In this script, I verify the validity of CSD association hits using the number of heterozygous females at each candidate.
# Cyril Matthey-Doret
# 04.12.2017


#=== LOADING ===#
# Loading packages
pack <- c("ggplot2","dplyr", "readr", "viridis", "reshape2", "tibble")
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
  filter( fisher > 3 & fisher == max(fisher))

# Subset CSD SNPs
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
for ( i in 1:nrow(csd_stat) ) {
  for (j in 1:nrow(csd_stat) ) {
    if ( i != j){
      csd_mat[i,j] <- length(csd_sum[csd_sum == 2^(i-1) + 2^(j-1) ])
    }else{
      csd_mat[i,j] <- length(csd_sum[csd_sum == 2^(i-1) ])
    }
  }
}
csd_plot <- melt(csd_mat)
ggplot(data=csd_plot,aes(Var1,Var2,fill=value)) + geom_tile() + scale_fill_viridis()

# Estimating mis-calling rate using portion of unlikely homo -> het transitions

snps <-sum_stat[,4:ncol(sum_stat)] %>% 
  # Excluding chromosome and bp infos
  filter(rowSums(.=='E')>0) %>%
  # Excluding SNPs that are fixed-homozygous
  t(.) %>%
  # Set names as rows and SNPs as columns
  data.frame(.) %>%
  rownames_to_column("Name") %>% 
  merge(indv,by="Name")

trgen <- data.frame(transition=rep(0,4),counts=rep(0,4))
i=1
for(tr in c("EE","EO","OO","OE")){
   Pgen <- substr(tr,1,1)  # Parent genotype
   Ogen <- substr(tr,2,2)  # Offspring genotype
   tmp <- snps %>%
      group_by(Family) %>%
      .[sapply(.[.$Generation == "F3",], 
               function(x) "F3" %in% x | Pgen %in% x)] %>%
      .[sapply(.[.$Generation == "F4",], 
               function(x) Ogen %in% x)] %>%
      ungroup() %>%
      select(starts_with('X')) %>%
      summarise(c = sum(. == Ogen))
   print(tmp)
   trgen[i,] <- c(tr,tmp[1,1])
   i <- i+1
}

test <- snps[,c(1:1000,356666, 356667, 356668)]
by(test[test$Generation=="F3",], INDICES = "Generation", function(x) 'E' %in% x)
