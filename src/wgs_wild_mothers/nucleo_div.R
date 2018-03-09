# Analysis of genome-wide nucleotidic diversity in wgs data from wild mothers
# Cyril Matthey-Doret
# 11.02.2018

packs <- c("ggplot2","dplyr", "gridExtra")
lapply(packs, require, character.only=T)

# Coordinates of CSD GWAS hits
csd_hits <- read.table("data/assoc_mapping/case_control/case_control_hits.tsv", header = T, stringsAsFactors = F)
# Filtering significant (p<10-3) hits. fischer column represents -log10 pval
csd_hits <- csd_hits %>%
  rename(CHROM=Chr) %>%
  filter(fisher >= 3)

# Anchored contigs containing CSD hits
csd_contig <- read.table("data/assoc_mapping/CSD_tig.tsv", stringsAsFactors = F)
colnames(csd_contig) <- c("CHROM", "START", "END")

# Nucleotidic diversity windows along wild mothers genomes
pi <- read.table("data/wgs_wild_mothers/stats/nucleo_div.windowed.pi", sep="\t", header=TRUE, stringsAsFactors = F)
pi <- pi[grep(pattern="chr", pi$CHROM),]

# CSD windows around significant SNPs
win_size <- 100000
csd_win <- csd_hits %>%
  select(CHROM, BP) %>%
  mutate(START=BP-win_size/2, END=BP+win_size/2, ID=row_number()) %>%
  select(-BP)


in_win <- function(rnum){
  tmp <- pi[pi$CHROM == csd_win[rnum, "CHROM"] & 
            pi$BIN_START > csd_win[rnum, "START"] &
            pi$BIN_END < csd_win[rnum, "END"],]
  if(nrow(tmp)){
    return(cbind(tmp,ID=csd_win$ID[rnum]))
  }
}

tmp <- sapply(1:nrow(csd_win), FUN=in_win, simplify = F)
csd_pi <- do.call(rbind,tmp)
csd_pi <- csd_pi[!duplicated(csd_pi),]

# Show Mb instead of basepairs
zoom=1000000

p0 <-ggplot(pi, aes(x=BIN_START/zoom, y=PI)) + 
  facet_grid(CHROM~., 
             space = 'free_x', 
             scales='free_x') + 
  geom_rect(data=csd_win, 
            inherit.aes = F, 
            aes(xmin=START/zoom, xmax=END/zoom, ymin=0, ymax=max(pi$PI)), alpha=0.3) +
  geom_point(col="blue", size=0.5) + geom_point(data=csd_pi, col="red", size=0.6) +
  xlab("Mb")
  
p1 <- ggplot(data=pi, aes(x=PI)) + 
  geom_histogram(fill="blue", aes(x=PI, ..density..), alpha=0.5, binwidth = 0.001) + 
  geom_histogram(data=csd_pi, fill="red", aes(x=PI, ..density..), alpha=0.3, binwidth = 0.001) +
  theme_bw()

p2 <- ggplot(data=pi, aes(x=1, y=PI)) + 
  geom_boxplot(fill="blue") + 
  geom_boxplot(data=csd_pi, fill="red", aes(x=2)) + 
  coord_flip() + theme_minimal()

grid.arrange(grobs=list(p0, p1,p2), layout_matrix=rbind(c(1,1,2),
                                                        c(1,1,2),
                                                        c(1,1,3)))
#pi <- pi[pi$CHROM=="chr3",]
#pi_cwt <- cwt(input = pi$PI, noctave = 4, nvoice = 8, plot = T)
#plot(pi_cwt[,1])


# Compute maximum PI in each CSD window to rule out those with only low values
pi$ID=as.factor(0);csd_pi$ID <- as.factor()
ggplot(data=pi, aes(x=ID, y=PI)) + 
  geom_boxplot(fill="blue") + 
  geom_boxplot(data=csd_pi, fill="red") + 
  coord_flip() + theme_minimal()
