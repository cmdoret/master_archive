# General statistics, per individual and per locus to have an overview of the data.
# Cyril Matthey-Doret
# 14.07.2017

# Per individual

indv <- read.table("../../data/ploidy/thresholds/fixed",header = T)

mother <- indv[indv$Generation=="F3",]
offspring <- indv[indv$Generation=="F4",]
prop_males <- function(mother){
  # Given a mother (row of dataframe), compute the proportion of males
  # among its dipoid offspring from the "offspring" dataframe and return it.
  diplo <- offspring[offspring$Ploidy=='D' & offspring$Family==mother["Family"],]
  male2N <- diplo[diplo$Sex=="M",]
  return(nrow(male2N)/nrow(diplo))
}
mother$prop_males = apply(mother,MARGIN = 1,FUN=prop_males)
# Proportion of males among dipoid offspring for each mother

hist(mother$prop_males,breaks=10, main="Proportion of males in diploid offspring per mother", xlab="Proportion of males")
hist(indv$N_SITES, breaks=20, main="Number of sites per individual (both 1n and 2n)", xlab="Number of sites")
hist(indv$MEAN_DEPTH,breaks=20, main="Mean depth per individual (both 1n and 2n)", xlab="Mean depth")

# Per locus

