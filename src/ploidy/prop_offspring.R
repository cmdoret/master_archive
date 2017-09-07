# This script is used to visualize offspring proportions per female.
# It takes an input list of individuals with sex, family, generation 
# and ploidy informations.
# Cyril Matthey-Doret
# 19.06.2017

load_packages <- function(){all(require(dplyr), require(tidyr), require(ggplot2))}
require(reshape2)
if(!suppressMessages(require(tidyverse))){
  stopifnot(load_packages())
}

in_path <- commandArgs(trailingOnly = T)[1]  # Taking ploidy list as CL argument
in_table <- read.table(in_path, header=T)

off_prop <- in_table %>%
  filter(Generation=='F4') %>%
  select(Family, Sex, Ploidy) %>%
  group_by(Family) %>%
  summarise(`Haploid male`=sum(Ploidy=='H'), 
            `Diploid male`=sum(Ploidy=='D' & Sex=='M'), 
            `Female`=sum(Sex=='F'))

off_prop <- melt(off_prop, id.vars = "Family", value.name = "Count", 
     measure.vars = c("Haploid male", "Diploid male", "Female"))

pdf(paste0('reports/lab_book/ploidy_per_fam/', basename(in_path), ".pdf"))
ggplot(data=off_prop, aes(x=Family, y=Count, fill=variable)) + geom_bar(stat='identity')
dev.off()
write.table(store_fam_ploidy, file = paste0('reports/lab_book/ploidy_per_fam/', basename(in_path)), row.names = F)
