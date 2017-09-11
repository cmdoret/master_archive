# This script uses the genomic output from STACK's populations module
# to identify potential individuals with terminal fusion automixis (as 
# opposed to central fusion automixis). This is done by measuring 
# heterozygosity rates in the centromere region inferred by chrom_types.R
# Cyril Matthey-Doret
# 29.08.2017


pack <- c("stringr","dplyr", "readr", "zoo", "ggplot2", "ggjoy", "reshape2")
lapply(pack, require, character.only = TRUE)

#==== LOADING DATA ====#

# Path to genotype file. Encoded as 3 letters: E (het), O (hom) and M (missing)
geno_path <- "../../data/assoc_mapping/grouped_geno_EOM.tsv"
gen <- read_tsv(geno_path, col_names = T)
# Loading centromere list
centro <- read_tsv("../../data/assoc_mapping/centro/centrolist.tsv")
# Loading list of individuals
indv <- read_tsv("../../data/individuals.tsv")
# Keeping only individuals present in the genotype file and excluding mothers
indv <- indv %>% 
  filter(Name %in% colnames(gen)[4:dim(gen)[2]]) %>%
  filter(Generation=='F4')

#==== PROCESSING ===#

# Excluding unordered contigs
gen <- gen %>% 
  filter(str_detect(string=Chr, pattern="chr.*"))

# Optional, only keeping large families for visualization
big_fam <- indv %>%
  group_by(Family) %>%
  filter(n()>12)
big_fam <- unique(big_fam$Family)

#==== ANALYSIS ====#

# Adding centromere positions to genotype file
gen <- gen %>% 
  group_by(Chr) %>%
  merge(centro, by="Chr")
# Using centromere position to compute distance from centromere
gen <- gen %>%
  mutate(centroD=abs(gen$pos - gen$BP)) %>%
  select( -pos, -val)
# Removing all sites which are missing for all offspring
gen <- gen %>% filter_at(vars(indv$Name), any_vars(str_detect(string=.,pattern="[OE]")))

#' slide_het
#' Computes the proportion of heterozygous sites in a sliding window.
#' @param seq A vector of characters, containing E (heterozygous) or O (homozygous) information.
#' @param dist A vector of distance to centromere associated to the sites in seq.
#' @param w The desired width of each window in number of sites.
#' @param step The step between each window, in number of sites.
#'
#' @return a zoo series with the proportion of heterozygous sites in each window with associated 
#' distance to centromere. The distance to centromere used is at the center of each window.
slide_het <- function(seq,dist,w=30,step=1){
  
  zoo_seq <- zoo(seq,order.by = dist)
  win_het <- rollapply(zoo_seq,width = w, by=step, align='center', function(x){sum(x=='E')/sum(x %in% c('E','O'))})
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
full_win$centro_dist <- as.numeric(full_win$centro_dist)

genofull= data.frame()
for(size in seq(50000,5000000,50000)){
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
# Merge chromosomes
genofull$Chr <- 'all_chr';full_win$Chr <- 'all_chr'
#==== VISUALISATION ====#
# Indexing individuals separately for each family
uniq_win <- unique(full_win[,c("Family","Name")])
# Creating smaller dataframe for indexing
uniq_win <- uniq_win %>% group_by(Family) %>% mutate(colindex=row_number())
# merging dataframes for index
full_win <- merge(full_win, uniq_win, by=c("Family","Name"))
# Computing normalized heterozygosity (not used atm)
full_win <- full_win %>% 
  group_by(Chr,Family) %>%
  mutate(norm_het = (Het. - mean(Het.))/sd(Het.))
# Computing median heterozygosity per family (could be used to classify CFA/TFA)
full_win <- full_win %>% group_by(Name) %>% mutate(medhet=median(Het., na.rm=T))
ggplot(data=full_win, aes(x=Name, y=Het., fill=as.factor(colindex))) + geom_boxplot() + 
  facet_wrap(~Family, scales='free_x') + guides(fill=F)
ggplot(data=full_win, aes(x=centro_dist, y=Het., col=as.factor(colindex),group=Name)) + 
  geom_line(stat='smooth', method='loess', se=F) + guides(col=FALSE) + facet_wrap(~Family)
# No correlation between number of loci in centromeric region and proportion of het.
smoothScatter(genofull$Num.Loci, genofull$Het.)
famorder <- cen_chr %>%
  group_by(Family) %>%
  summarise(avg=mean(Het.))
cen_chr$Family <- factor(cen_chr$Family, ordered = T, levels = famorder$Family[order(famorder$avg)])
ggplot(data=genofull, aes(x=Het., y=Family)) + geom_joy()
genofull <- genofull %>% 
  group_by(Chr, Family) %>%
  mutate(norm_het=(Het.-mean(Het., na.rm=T))/sd(Het., na.rm=T))

ggplot(data=genofull[genofull$Family %in% big_fam,], aes(x=centrosize, y=Het.,col=Name, weight=Num.Loci)) + stat_smooth(method='loess') + 
  facet_grid(Family~Chr)+ guides(col=FALSE)
