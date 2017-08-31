# This script uses the genomic output from STACK's populations module
# to identify potential individuals with terminal fusion automixis (as 
# opposed to central fusion automixis). This is done by measuring 
# heterozygosity rates in the centromere region inferred by chrom_types.R
# Cyril Matthey-Doret
# 29.08.2017


pack <- c("stringr","dplyr", "readr")
lapply(pack, require, character.only = TRUE)
centrosize <- 2000000

#==== LOADING DATA ====#

# Path to folder containing STACKS populations output files
pop_path <- "../../data/populations/grouped_d-3_r-80"
# File with sample names ordered as in genomic output
sampleID <- read_tsv(paste(pop_path, 'batch_0.sumstats.tsv', sep='/'),
                  n_max = 2, col_names = F)
# Genomic output from populations
gen <- read_tsv(paste(pop_path, "batch_0.genomic.tsv", sep='/'), 
                skip = 1, col_names = F)
# Loading centromere list
centro <- read_tsv("../../data/assoc_mapping/centro/centrolist.tsv")

#==== PROCESSING ===#

# Concatenating sample names into a single vector
sampleID <- append(unlist(str_split(sampleID[1,2], ',')),
                   unlist(str_split(sampleID[2,2], ',')))
# Meaningful column names
colnames(gen)[1:3] <- c("Locus.ID","Chr","BP")
gen <- gen %>%
  filter(str_detect(string= Chr,pattern = "chr.*"))
# Replacing numeric headers with sample names
colnames(gen)[4:length(colnames(gen))] <- sampleID
# If populations error yielded high number, replace with missing code
gen[,sampleID] <- apply(gen[,sampleID], 2, function(x) ifelse(x>10,yes = 0,no=x))
# Get sample columns into character format
gen[,sampleID] <- apply(gen[,sampleID], 2, as.character)
# Constructing genotype dictionary to match numeric genotype with hetero/homozygous/missing
geno <- list('0'='M')
for(g in as.character(1:10)){geno[[g]] <- ifelse(g %in% c('1','5','8','10'), yes='O', no='E')}
# Replacing numeric encoding with matching genotype letter
vec_geno <- Vectorize(function(x) geno[[as.character(x)]])
gen[,sampleID] <- lapply(gen[,sampleID], FUN=vec_geno)
#gen <- gen %>% filter_at(vars(sampleID), any_vars(.=='E'))

centro$start <- centro$pos - centrosize
centro$end <- centro$pos + centrosize

#==== ANALYSIS ====#

for ( chrom in unique(gen$Chr)){
  chr_start <- centro$start[centro$Chr == chrom]
  chr_end <- centro$end[centro$Chr == chrom]
  cen_chr <- gen %>% 
    filter(Chr == chrom & BP > chr_start & BP < chr_end) %>%
    summarise_at(vars(sampleID), function(x){sum(x=='E')/sum(x %in% c('E','O'))})
}

#==== VISUALISATION ====#

geno[[as.character(gen[1,'C209'])]]
