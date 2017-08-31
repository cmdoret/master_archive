# This script uses the genomic output from STACK's populations module
# to identify potential individuals with terminal fusion automixis (as 
# opposed to central fusion automixis). This is done by measuring 
# heterozygosity rates in the centromere region inferred by chrom_types.R
# Cyril Matthey-Doret
# 29.08.2017


pack <- c("stringr","dplyr", "readr", "zoo", "ggplot2", "ggjoy", "reshape2")
lapply(pack, require, character.only = TRUE)
centrosize <- 1000000

#==== LOADING DATA ====#

# Path to folder containing STACKS populations output files
geno_path <- "../../data/assoc_mapping/grouped_geno_EOM.tsv"
# Genomic output from populations
gen <- read_tsv(geno_path, col_names = T)
# Loading centromere list
centro <- read_tsv("../../data/assoc_mapping/centro/centrolist.tsv")
# Loading list of individuals
indv <- read_tsv("../../data/individuals")
indv <- indv %>% 
  filter(Name %in% colnames(gen)[4:dim(gen)[2]]) %>%
  filter(Generation=='F4')
#==== PROCESSING ===#

gen <- gen %>% 
  filter(str_detect(string=Chr, pattern="chr.*"))
centro$start <- centro$pos - centrosize
centro$end <- centro$pos + centrosize

#==== ANALYSIS ====#

gen <- gen %>% 
  group_by(Chr) %>%
  merge(centro, by="Chr")
gen <- gen %>%
  mutate(centroD=abs(gen$pos - gen$BP)) %>%
  select(-start, -end, -pos, -val)

# !!! CHECK THIS LINE !!!
gen <- gen %>% filter_at(vars(indv$Name), any_vars(str_detect(string=.,pattern="[OE]")))

slide_het <- function(seq,dist,w=50,step=1){
  zoo_seq <- zoo(seq,order.by = dist)
  win_het <- rollapply(zoo_seq,width = w, by=step, function(x){sum(x=='E')/sum(x %in% c('E','O'))})
  return(win_het)
}

full_win <- data.frame()
for(chrom in unique(gen$Chr)){
  chr_gen <- gen[gen$Chr==chrom,]
  win_chr <- data.frame(lapply(chr_gen[,indv$Name],  function(x) slide_het(x, chr_gen$centroD)))
  if(nrow(win_chr)>0){
    win_chr$Chr <- chrom
    win_chr$centro_dist <- rownames(win_chr);rownames(win_chr) <- NULL
    full_win <- rbind(full_win, win_chr)
  }
}
full_win <- melt(full_win, measure.vars = indv$Name, id.vars = c("Chr","centro_dist"), variable.name="Name",value.name="Het.")
full_win <- merge(full_win, indv, by="Name")

genofull= data.frame()
for(size in c(100000,3000000)){
  for ( chrom in unique(gen$Chr)){
    chr_start <- centro$pos[centro$Chr==chrom] - size
    chr_end <- centro$pos[centro$Chr==chrom] + size
    cen_chr <- gen %>% 
      filter(Chr == chrom & BP > chr_start & BP < chr_end) %>%
      select(indv$Name[indv$Generation=='F4'])
    chr_prop <- cen_chr %>% 
      summarise_all(function(x){sum(x=='E')/sum(x %in% c('E','O'))})
    chr_loci <- cen_chr %>% 
      summarise_all(function(x){sum(x %in% c('E','O'))})
    cen_chr <- data.frame(t(rbind(chr_prop, chr_loci)))
    colnames(cen_chr) <- c("Het.", "Num.Loci")
    cen_chr <- cen_chr %>%
      tibble::rownames_to_column(var = "Name")%>%
      merge(indv, by="Name")
    cen_chr <- cen_chr %>%
      mutate(centrosize=size, Chr=chrom)
    genofull <- rbind(genofull, cen_chr)
  }
  print(size)
}
#==== VISUALISATION ====#

ggplot(data=full_win, aes(x=centro_dist, y=Het., group=Name, col=Family)) + geom_line(stat='smooth', method='loess',alpha=0.3, se=F) + 
  facet_grid(~Chr) + guides(col=FALSE)

# No correlation between number of loci in centromeric region and proportion of het.
smoothScatter(genofull$Num.Loci, genofull$Het.)
famorder <- cen_chr %>%
  group_by(Family) %>%
  summarise(avg=mean(Het.))
cen_chr$Family <- factor(cen_chr$Family, ordered = T, levels = famorder$Family[order(famorder$avg)])
ggplot(data=cen_chr, aes(x=Het., y=Family)) + geom_joy()
genofull <- genofull %>% 
  group_by(Chr, Family) %>%
  mutate(norm_het=(Het.-mean(Het., na.rm=T))/sd(Het., na.rm=T))

ggplot(data=genofull, aes(x=centrosize, y=Het., col=Family, group=Name, weight=Num.Loci)) + geom_line() + 
  facet_grid(.~Chr)+ guides(col=FALSE)
