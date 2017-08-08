hom_filt <- function(fam){
  # Function for filtering loci homozygous in mother of given family
  mother <- as.character(indv$Name[indv$Family==fam & indv$Generation=='F3'])
  # Using mother's genotype from snps file
  snp_gz <- gzfile(paste0('../../data/sstacks/',fam,'/',mother,'.snps.tsv.gz'),"rt")
  # Uncompresing file (lengthy process)
  snp_mother <- read.table(snp_gz); close(snp_gz)
  # keeping all loci where all SNPs are hom
  hom_loc <- snp_mother %>% group_by(V3) %>% filter(all(V5=='O'))  
  # Find loci which have only homozygous positions in pstacks snps file
  hom_ID_local <- unique(hom_loc$V3)
  match_gz <- gzfile(paste0('../../data/sstacks/',fam,'/',mother,'.matches.tsv.gz'),"rt")
  # uncompressing...
  match_mother <- read.table(match_gz); close(match_gz)
  hom_ID_global <- unique(match_mother$V3[match_mother$V4 %in% hom_ID_local])
  # Finding catalog locus IDs matching sample locus IDs using sstacks matches file
  return(hom_ID_global)
}

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
# Customized operator for convenience

# Genome statistics
# phi_path <- commandArgs(TrailingOnly=T)[1]

if(FALSE){
  rec_hom <- data.frame()
  sum_stat <- data.frame()
  for(fam in list.dirs(stat_path)[2:length(list.dirs(stat_path))]){  # Excluding first dir (parent)
    tmp_stat <- read.csv(paste0(fam,'/batch_0.sumstats.tsv'),header=T,skip=2,sep='\t')
    tmp_stat$fam <- rep(basename(fam))
    try(mother_hom <- hom_filt(basename(fam)))
    rec_hom <- rbind(rec_hom, tmp_stat[tmp_stat$Locus.ID %in% mother_hom,c("Chr","BP","fam")])
    tmp_stat <- tmp_stat[tmp_stat$Locus.ID %not in% mother_hom,]
    sum_stat <- rbind(sum_stat, tmp_stat)
    mother_hom <- NULL
  }
}

if(FALSE){
  rec_hom <- read.csv("../../data/SNP_lists/fixed_hom_mother.txt",header=T)
  sum_stat <- data.frame()
  for(fam in list.dirs(stat_path)[2:length(list.dirs(stat_path))]){  # Excluding first dir (parent)
    tmp_stat <- read.csv(paste0(fam,'/batch_0.sumstats.tsv'),header=T,skip=2,sep='\t')
    tmp_stat$fam <- rep(basename(fam))
    sum_stat <- rbind(sum_stat, tmp_stat)
    mother_hom <- NULL
  }
  
  rec_hom$EXTRACOL <- 1 # use a column name that is not present among 
  # the original data.frame columns
  merged <- merge(sum_stat,rec_hom,all=TRUE, by.x=c("Chr","BP","fam"), by.y=c("contig","bp","family"))
  
  sum_stat <- merged[is.na(merged$EXTRACOL),]
  sum_stat$EXTRACOL <- NULL
}

