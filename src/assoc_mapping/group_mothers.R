# This script groups mothers based on the proportion of males among their diploid offspring.
# These different groups will have different combination of heterozygous CSD loci.
# Cyril Matthey-Doret
# 30.07.2017

library(dplyr)
indv <- read.table("../../data/ploidy/thresholds/fixed", header=T)

diplo <- indv[indv$Generation=="F4" & indv$Ploidy=='D',]
off_comp <- data.frame(table(diplo$Family, diplo$Sex))
colnames(off_comp) <- c("Family", "Sex", "Count")

off_comp <- off_comp %>%
  group_by(Family) %>%
  mutate(prop_males=round(Count[Sex=='M']/sum(Count),3)) %>%
  filter(Sex=='M') %>%
  select(-Sex)

hist(off_comp$prop_males,breaks=10, main="Male proportion in each family", xlab="Proportion of males", ylab="Number of families")
print(off_comp)
plot(off_comp$prop_males,off_comp$Count, xlab="Proportion of males", 
     ylab="Number of 2N offspring", main="Proportion of males versus total diploid offpsring")

# Arbitrary cutoff needed here