# Trying to summarise the importance of different CSD loci with tables/viualisations
# Cyril Matthey-Doret
# 14.12.2017


#=== LOADING ===#
# Loading packages
pack <- c("ggplot2","dplyr", "readr", "viridis", "reshape2", "tibble", "eulerr", "VennDiagram", "gridExtra")
if (all(unlist(lapply(pack, require, character.only = TRUE)))){
  print("Packages loaded !")
}
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

# Subset CSD SNPs and flag het females
fem_het<- sum_stat %>%
  right_join(., best, by=c("Chr","BP")) %>%
  select(2,3,one_of(indv$Name[indv$Sex=='F'])) %>%
  mutate_at(funs(ifelse(.=='E',1,0)),.vars = vars(-Chr,-BP))

# Subset CSD SNPs and flag het males
male_het <- sum_stat %>%
  right_join(., best, by=c("Chr","BP")) %>%
  select(2,3,one_of(indv$Name[indv$Sex=='M'])) %>%
  mutate_at(funs(ifelse(.=='E',1,0)),.vars = vars(-Chr,-BP))


# Number of loci
L <- nrow(fem_hom)

comb_csd <- list()
tot_comb <- 0
for (i in 1:L) {
  comb_csd[[i]] <- combn(c(1:L), i)
  tot_comb <- tot_comb + ncol(comb_csd[[i]])
}
euler_dat <- list(name = rep("", tot_comb), 
                  fem  = rep(0, tot_comb), 
                  male = rep(0, tot_comb))
#venn_names <- rep("",tot_comb)
iter <- 1
for (n in comb_csd){
  for (i in 1:ncol(n)) {
    incl_loci <- n[, i]
    euler_dat$fem[iter] <- sum(colSums(fem_het[incl_loci, 3:dim(fem_het)[2]]) == 0)
    euler_dat$male[iter] <- sum(colSums(male_het[incl_loci, 3:dim(male_het)[2]]) > 0)
    euler_dat$name[iter] <- paste(LETTERS[incl_loci], sep="&",collapse = "&")
    prefix <- ifelse(length(incl_loci) > 1, yes="n", no="area")
    #venn_names[iter] <- paste0(prefix, paste(incl_loci, collapse=""))
    iter <- iter + 1
  }
}

# Proprtional areas (Euler)
names(euler_dat$fem) = names(euler_dat$male)  <- euler_dat$name
p_fem <- plot(euler(euler_dat$fem,shape = "ellipse"),quantities = T,fill_alpha = 0.2, 
     main="Females hom. at all CSD peak")
p_male <- plot(euler(euler_dat$male,shape = "ellipse"),quantities = T,fill_alpha = 0.2, 
              main="Males het. at any CSD peak")

grid.arrange(p_fem,p_male,nrow = 1)

# Fixed areas (Venn)
#for (i in 1:L){
#  euler_store[i] <- sum(euler_store[grepl(venn_names, pattern = paste0(".*",i,".*"))])
#names(euler_store) <- venn_names
#}
#do.call(draw.quad.venn,as.list(euler_store))

rowSums(male_csd[,3:dim(male_csd)[2]] == 1)
