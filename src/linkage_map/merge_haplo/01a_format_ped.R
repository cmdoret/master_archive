# Selecting families that will be used to build the linkage map.
# Small families and those where the mother is missing are excluded.
# Cyril Matthey-Doret
# 22.05.2018

library(dplyr);library(readr);library(stringr)

args <- commandArgs(trailingOnly = TRUE)

# Tabular file containing sample informations (ploidy file)
samples <- read_tsv(args[1], col_names=T, col_types=cols())

# Only including families with at least SIZE individuals 
# Before diploidization and including mother
SIZE <- as.numeric(args[2])

# Matrix to write LEPMAP3-compatible pedigree file
ped <- matrix(ncol=(nrow(samples)+2), nrow=6)

# All lines must start with CHROM POS
ped[, 1] <- "CHROM"
ped[, 2] <- "POS"

# First line contains family names
ped[1, 3:ncol(ped)] <- samples$Family

# Second line contains sample ID
ped[2, 3:ncol(ped)] <- samples$Name

# Third line contains father ID
# Fathers set to 0
ped[3:4, (which(samples$Generation=='F3' & samples$Sex=='M') + 2)] <- 0
# Fill indices of F4 with their respective (synthetic) father ID
for (fam in unique(samples$Family)){
  ped[3, (which(samples$Generation=='F4' & 
                  samples$Family==fam) + 2)] <- samples$Name[samples$Generation=='F3' &
                                                               samples$Family==fam &
                                                               samples$Sex=='M']
}

# Fourth line contains mother ID
# Mothers set to 0
ped[3:4, (which(samples$Generation=='F3' & samples$Sex=='F') + 2)] <- 0
# Fill indices of F4 with their respective mother ID
for (fam in unique(samples$Family)){
  ped[4, (which(samples$Generation=='F4' & 
                samples$Family==fam) + 2)] <- samples$Name[samples$Generation=='F3' &
                                                           samples$Family==fam &
                                                           samples$Sex=='F']
}

# Fifth line contains sex information
ped[5, 3:ncol(ped)] <- recode(samples$Sex,M=1, F=2)
    
# Sixth line could contain a phenotype (here: placeholder 0's)
ped[6, 3:ncol(ped)] <- 0

# Order and filter samples as in VCF, if provided
if(!is.na(args[4])){
  # If samples in the VCF file are specified, get list of names to filter/order pedigree samples
  con_read <- pipe(paste0("grep '^#[^#]' ", args[4], " | sed 's/.*FORMAT.//'"))
  lep_samples <- scan(con_read, what='character')
  close(con_read)
  colnames(ped) <- ped[2,]
  ordered_cols <- append(c("CHROM", "POS"), lep_samples)
  ped <- ped[,order(match(colnames(ped), ordered_cols, nomatch = NA), na.last=NA)]
}

write.table(ped, args[3], col.names = F, row.names = F, quote = F, sep='\t')


# The likelihoods will be provided on a separate file given as parameter posteriorFile or vcfFile
