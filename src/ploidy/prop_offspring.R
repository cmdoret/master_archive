# This script is used to visualize offspring proportions per female.
# It takes an input list of individuals with sex, family, generation 
# and ploidy informations.
# Cyril Matthey-Doret
# 19.06.2017

require(dplyr); require(tidyr); require(ggplot2); require(reshape2)

in_path <- commandArgs(trailingOnly = T)[1]  # Taking ploidy list as CL argument
in_table <- read.table(in_path, header=T)

off_prop <- in_table %>%
  filter(Generation=='F4') %>%
  select(Family, Sex, Ploidy) %>%
  group_by(Family) %>%
  summarise(`Haploid son`=sum(Ploidy=='H'), 
            `Diploid son`=sum(Ploidy=='D' & Sex=='M'), 
            `Daughter`=sum(Sex=='F'))

off_prop_long <- melt(off_prop, id.vars = "Family", value.name = "Count", 
     measure.vars = c("Haploid son", "Diploid son", "Daughter"))
off_prop_long$Family <- as.character(off_prop_long$Family)
sum_count <- off_prop_long %>%
  group_by(Family) %>% 
  summarise(tot=sum(Count))
off_prop_long$Family <- factor(off_prop_long$Family, levels=sum_count$Family[order(sum_count$tot)], ordered=T)

pdf(paste0('reports/lab_book/ploidy_per_fam/', "off_prop.pdf"))
ggplot(data=off_prop_long, aes(x=Family, y=Count, fill=variable)) + geom_bar(stat='identity')
dev.off()
write.table(off_prop, file = paste0('reports/lab_book/ploidy_per_fam/', 'off_prop.tsv'), sep='\t', quote=F, 
                                    row.names = F)
