# Reconstructs asexual mothers' genotype from their offspring's 
# CMD
# 06.06.2018

#### GENOTYPE ENCODING ####

# As of STACKS 1.48, genotypes are encoded as follows
# In the populations "genomic" output file
#    A  C  G  T
# A  1  2  3  4
# C  2  5  6  7
# G  3  6  8  9
# T  4  7  9  10

# Hence, homo/hemizygous genotypes are 
# A:1 C:5 T:8 G:10

# Packages
packs <- c("dplyr","readr","ggplot2", "tidyr", "stringr", "optparse")
packs <- sapply(packs, function(x) 
  suppressPackageStartupMessages(library(x, quietly=T, character.only=T)))

#### PARSE CL ARGUMENTS ####

option_list = list(
  make_option(c("-o", "--out"), 
              type="character", 
              default=NA, 
              help="Output file where to store reconstructed mothers genotypes.", 
              metavar="character"),
  make_option(c("-f", "--fam"), 
              type="character", 
              help="The 'family' file. Contains Names, Families, Sex and Mother ID."),
  make_option(c("-l","--plo"), 
              type="character", 
              default=NA, 
              help="The ploidy file. Contains Names, Families, Sex and Ploidies."),
  make_option(c("-p","--pop"), 
              type="character", 
              help="The stacks populations output folder. The 'genomic' and 'sumstats' files will be used and should be prefixed by batch_1."),
  make_option(c("-m","--mfreq"), 
              type="integer", 
              default=10, 
              help="Minimum proportion (integer, in percentage) of offspring in which an allele must be called to be considered. Default: 10"),
  make_option(c("-L","--LEPMAP"), 
              action="store_true",
              default=F,
              help="If option is specified, the output will be stored as LEPMAP compatible format with a pedigree header")

)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# check all paths are provided
if (any(is.na(opt$out), is.na(opt$pop), is.na(opt$fam), is.na(opt$plo))) {
 stop("You need to provide paths for all input and output files. See script usage (--help)")
}


# Minimum allele frequency required to include in mother
min_freq = opt$mfreq / 100
if (min_freq > 1 | min_freq < 0) {
  stop("You need to provide a minimum allele frequency between 0 and 100 (-m | --mfreq). See script usage (--help)")
}

#### LOAD ####

# Family: Used to retrieve mothers' ID
fam <- read_tsv(
    pipe(paste0("cut -f1-4 ", opt$fam)), 
    col_names = T, col_types='cccc') %>% 
  select(Family, Parent_id) %>%
  distinct(Family, Parent_id)

# Read ploidy file
ploid <- read_tsv(opt$plo, col_names = T, col_types = cols()) %>%
  # Only keep relevant variables
  select(Name, Sex, Family, Generation, Ploidy)

# Read sample names into a vector, in the same order as genotype matrix
clean_pop <- pipe(paste0(
  "grep '^# [MF]' ", opt$pop, "/batch_1.sumstats.tsv \\
  | sed 's/^# [MF]//' \\
  | tr '\n' ',' \\
  | tr -d '\t' \\
  | sed 's/,$//'"),open='r')

s_names <- scan(clean_pop, sep=',', what="character")

close(clean_pop)

# Geno matrix: Genotypes from all samples in analysis
geno <- read_tsv(paste0(opt$pop, "/batch_1.genomic.tsv"), 
                 col_names=F, skip=1, col_types=cols())

#### FORMAT ####

# Extract names of haploid sons
haplo_sons <- ploid$Name[ploid$Ploidy == 'H' & 
                         ploid$Generation == 'F4']
# And all diploid offspring
diplo_off <- ploid$Name[ploid$Ploidy == 'D' & 
                           ploid$Generation == 'F4']

# Name columns of genotype matrix after samples
colnames(geno) <- append(c("ID", "Chr", "BP"), s_names)

geno <- geno %>%
  # Reshape into "long" format
  gather(Name, Nuc, -Chr, -BP, -ID)

haplo <- geno %>% 
  # Only keep haploid samples
  filter(Name %in%  haplo_sons) %>%
  # Include Family information
  left_join(ploid %>% select(Name, Family), by="Name") %>%
  # Recode genotypes into nucleotides. missing and het -> N
  mutate(Nuc = case_when(
    Nuc == "1"  ~ "A",
    Nuc == "5"  ~ "C", 
    Nuc == "8"  ~ "G", 
    Nuc == "10" ~ "T", 
    TRUE ~ "N"))

diplo <- geno %>%
  # Only keep haploid samples
  filter(Name %in%  diplo_off) %>%
  # Include Family information
  left_join(ploid %>% select(Name, Family), by="Name")
  
#    A  C  G  T
# A  1  2  3  4
# C  2  5  6  7
# G  3  6  8  9
# T  4  7  9  10
# Vectors are 1-indexed in R -> shift genotypes IDs by 1 
#           1   2   3   4   5   6   7   8   9   10  11
geno1 <- c("N","A","C","G","T","C","G","T","G","T","T")
geno2 <- c("N","A","A","A","A","C","C","C","G","G","T")

# Recode genotypes into nucleotides in both haplotypes.
diplo1 <- diplo %>% mutate(Nuc = geno1[Nuc+1])
diplo2 <- diplo %>%
  mutate(Name = str_c(Name, "_2")) %>%
  mutate(Nuc = geno2[Nuc+1])

# Combine haploid and diploid samples
geno <- haplo %>% 
  bind_rows(diplo1) %>%
  bind_rows(diplo2)

#### RECONSTRUCT ####
# Reconstruct asexual mothers genotype using their offspring

geno <- geno %>%
  # Count occurrences of each genotypes per family
  group_by(ID, Chr, BP, Family, Nuc) %>%
  summarise(Fnuc = n()) %>%
  # Transform into proportions
  group_by(ID, Chr, BP, Family) %>%
  mutate(Fnuc = Fnuc / sum(Fnuc))

mother_g <- geno %>%
  # Group genotypes by family and SNP
  group_by(ID, Chr, BP, Family) %>% 
  # Sort by allele frequency
  arrange(desc(Fnuc)) %>%
  # Keep 2 most frequent alleles
  top_n(2, wt = Fnuc) %>% 
  # Mother gets assigned these alleles if above minimum frequency
  summarise(Nuc1  = ifelse(Fnuc[1] >= min_freq, 
              yes = Nuc[1], 
              no  = "N"), 
            Nuc2  = case_when(
              n() == 1 ~ Nuc1,
              Fnuc[2] >= min_freq ~ Nuc[2],
              Fnuc[2] <  min_freq ~ Nuc1)) %>%
  # Match families with desired mother ID
  ungroup() %>%
  left_join(fam, by="Family") 

#### FORMAT ####
# Reformat mother genotypes into wide format (better for linkage mapping/QTL softwares)

mother_out <- mother_g %>% 
  # Mother state: heterozygous (H) or homozygous (A/C/T/G)
  mutate(State = case_when(
    Nuc1 == "N"  ~ "-",
    Nuc1 == Nuc2 ~ Nuc1, 
    Nuc1 != Nuc2 ~ "H")) %>%
  # Spread into wide format (mothers as columns, SNPs as rows)
  select(-Family, -Nuc1, -Nuc2) %>%
  spread(Parent_id, State)
 

#### LEPMAP3 COMPATIBILITY ####

# Desired format:
# header: 6 lines pedigree
# body: SNP as rows, samples as tab-sep column.
# each genotype call = 10 space-separated 0<=int/float<=1 representing likelihoods
# 10 nums represent: AA AC AG AT CC CG CT GG GT TT
# For N samples, N+2 column. The 2 first columns being CHROM and POS
# genotypes as sum of nucleotides
# e.g. third genotype in LEPMAP order is AG = A + G = 1 + 4 = 5
# Therefore genusum[3] = 5
gensum=c(2, 3, 5, 9, 4, 6, 10, 8, 12, 16)
# Need array to reverse access genotype position from sum:
# e.g. sum2gen[5] = 3 -> AG
# Then transorm into binary to format genotype
# e.g. AG is third so format should be 0010000000 
# This can be done by taking 2 to the power of the 0-indexed position and mirroring the binary
# AG = 3 -> 3-1=2 -> bin((2^2)) = 0000000100 -> 0010000000

sum2bin=rep("X", max(gensum))
for ( i in 1:length(gensum)){
  # get binary value to represent genotype position
  # note: no need to mirror, first bit is already on the left
  # Format with spaces -> ready for output file
  sum2bin[gensum[i]] = paste(as.numeric(intToBits(2^(i-1))[1:length(gensum)]), sep=" ", collapse=" ")
}

if(opt$LEPMAP == T){
  mother_out <- mother_g %>% 
    # Encoding nucleotides as binary values -> sum = unique genotype
    # Note: any genotype >=20 will be N
    mutate(Nuc1 = recode(Nuc1, "A"=1, "C"=2, "G"=4, "T"=8, "N"=20)) %>%
    mutate(Nuc2 = recode(Nuc2, "A"=1, "C"=2, "G"=4, "T"=8, "N"=20)) %>%
    mutate(State = Nuc1 + Nuc2) %>%
    # Make genotype LEPMAP compatible
    mutate(State = ifelse( State < 20, 
                          yes = sum2bin[State], 
                          no = "0 0 0 0 0 0 0 0 0 0")) %>%
    # Spread into wide format (mothers as columns, SNPs as rows)
    select(-Family, -Nuc1, -Nuc2) %>%
    spread(Parent_id, State)
  
  mother_info <- mother_g %>%
    select(-ID, -Chr, -BP, Nuc1, Nuc2) %>%
    distinct(Family, Parent_id) %>%
    mutate(Sex="F")
  father_info <- mother_info %>%
    mutate(Parent_id = paste0(Parent_id,  "_DUP", sep="")) %>%
    mutate(Sex="M")

  parent_info <- mother_info %>%
    bind_rows(father_info) %>% 
    arrange(Parent_id)
  
  # Building pedigree header
  ped <- as.data.frame(matrix(nrow=6, ncol=(nrow(parent_info)+2)))
  ped[,1] <- "CHROM"
  ped[,2] <- "POS"
  ped[1, 3:ncol(ped)] <- parent_info$Family
  ped[2, 3:ncol(ped)] <- parent_info$Parent_id
  # Father, Mother and placeholder for phenotype all set to 0
  ped[c(3, 4, 6), 3:ncol(ped)] <- 0
  # Sex encoded as int
  ped[5, 3:ncol(ped)] <- recode(parent_info$Sex, M=1, F=2)
  colnames(ped) <- ped[2,]
  colnames(ped)[1:2] <- c("Chr", "BP")
  
  # Duplicating genotype calls for synthetic fathers
  father_body <- mother_out %>% 
    select(-ID, -Chr, -BP) %>% 
    rename_all(funs(. = paste0(., "_DUP")))
  
  # Combine all parents genotypes into asingle table
  geno_lm3 <- mother_out %>%
    select(-ID) %>%
    bind_cols(father_body) %>%
    # Order columns in the same order as ped
    select(order(match(names(.),names(ped)))) %>%
    mutate(BP = as.character(BP))

  mother_out <- ped %>% bind_rows( geno_lm3)
}

#### WRITE ####

write_tsv(mother_out, opt$out )
