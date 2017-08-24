# This script groups mothers based on the proportion of males among their diploid offspring.
# These different groups will have different combination of heterozygous CSD loci.
# Cyril Matthey-Doret
# 30.07.2017

#==== LOADING DATA ====#
library(dplyr);library(ggplot2)
indv <- read.table("../../data/ploidy/thresholds/fixed", header=T)
scenario <- 2  # Number of CSD loci considered

#==== PROCESSING ====#
diplo <- indv[indv$Generation=="F4" & indv$Ploidy=='D',]
off_comp <- data.frame(table(diplo$Family, diplo$Sex))
colnames(off_comp) <- c("Family", "Sex", "Count")

off_comp <- off_comp %>%
  group_by(Family) %>%
  mutate(prop_males=round(Count[Sex=='M']/sum(Count),3)) %>%
  filter(Sex=='M') %>%
  select(-Sex)

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

km_output <- kmeans(off_comp$prop_males,centers = ngroups)
off_comp$cluster <- km_output$cluster
#ggplot(data=off_comp, aes(x=prop_males, fill=as.factor(cluster))) + geom_histogram() +
#  ggtitle("Male proportion per family") + xlab("Proportion of males") + theme_bw() +
#  scale_fill_discrete(name="Heterozygous loci", labels=c("a","a+b", "b")) + ylab("Number of families") +
#  ylab("Number of families") + geom_vline(data=data.frame(km_output$centers),
#                                          aes(xintercept=km_output.centers), col='red', lty=2)


plot(off_comp$prop_males,off_comp$Count, xlab="Proportion of males", 
     ylab="Number of 2N offspring", main="Proportion of males versus total diploid offpsring")

# Arbitrary cutoff needed here