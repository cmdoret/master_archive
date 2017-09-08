# This script groups mothers based on the proportion of males among their diploid offspring.
# These different groups will have different combination of heterozygous CSD loci.
# Number of diploid offspring is inferred by extrapolating the haploidy rate among sequenced 
# individuals to all (non-sequenced) offspring in the family.
# Cyril Matthey-Doret
# 30.07.2017

#==== LOADING DATA ====#

library(dplyr);library(ggplot2)
args_list <- commandArgs(trailingOnly = T)
scenario <- as.numeric(args_list[1])  # Number of CSD loci considered
indv <- read.table(args_list[2], header=T)  # Path to input file
out <- args_list[3]  # Output path

# Loading additional files
fam_list <- read.table(pipe("cut -f2,4 ../../data/families.txt"), header=T)
# Total number of non-sequenced offspring in each family
tot_off <- read.table("../../data/total_offspring.tsv", header=T)

#==== PROCESSING ====#
fam_list <- unique(fam_list)  # Used for mother ID - Family name correspondance
tot_off <- merge(tot_off, fam_list, by.x='mother', by.y='parent.id', all=F)
male_seq <- indv[indv$Sex=='M',]  # Sequenced offspring
male_pl <- data.frame(table(male_seq$Family, male_seq$Ploidy))
colnames(male_pl) <- c("Family", "Ploidy", "Count")

male_pl <- male_pl %>%
  group_by(Family) %>%
  mutate(prop_2n=round(Count[Ploidy=='D']/sum(Count),3), 
          tot_M=sum(Count)) %>%
  filter(Ploidy=='D') %>%
  select(-Ploidy, -Count)

diplo_off <- tot_off %>%
  merge(male_pl, by="Family",all=F) %>%
  mutate(infer_2n=sons*prop_2n) %>%
  select(-unknown, -prop_2n, -tot_M, -sons) %>%
  mutate(prop_male=round(infer_2n/(daughters+infer_2n),3))

#=== CLUSTERING ===#
# Expected number of categories is a function of the number of CSD loci
ngroups <- 2^scenario-1
# exp_prop <- 1/(2^scenario)

# Using Pascal's triangle to get count of genotypes 
# with n het loci where 0 < n <= scenario
# Note: n < 0 since 0 het loci lead to male
#x <- 1
#for (i in 1:scenario) { x <- c(0, x) + c(x, 0)}

# Number of centers for k-means clustering
#centers <- sum(x[2:length(x)])  # First value is only homozygous loci

km_output <- kmeans(diplo_off$prop_male,centers = ngroups)
diplo_off$cluster <- km_output$cluster
ggplot(data=diplo_off, aes(x=prop_male, fill=as.factor(cluster))) + geom_histogram() +
 ggtitle("Male proportion per family") + xlab("Proportion of males") + theme_bw() +
 scale_fill_discrete(name="Heterozygous loci", labels=c("a","a+b", "b")) + ylab("Number of families") +
 ylab("Number of families") + geom_vline(data=data.frame(km_output$centers),
                                         aes(xintercept=km_output.centers), col='red', lty=2)

write.table(off_comp, out, sep='\t', row.names=F, quote=F)
