# This script is used to differentiate between telocentric and metacentric
# chromosomes. This is done by using the density of heterozygous sites in the
# mother that become homozygous in offspring along the chromosome.
# 30.07.2017
# Cyril Matthey-Doret

library(ggplot2); library(dplyr);library(viridis)


# Genome statistics
stat_path <- '../../data/populations/d-3_r-80/'
indv <- read.table('../../data/individuals', header=T)
grouped <- "T"
# phi_path <- commandArgs(TrailingOnly=T)[1]

# removing SNPs homozygous in all females (consequently in mothers)
if(grouped=='F'){
  sum_stat <- data.frame()
  for(fam in list.dirs(stat_path)[2:length(list.dirs(stat_path))]){  # Excluding first dir (parent)
    tmp_stat <- read.csv(paste0(fam,'/batch_0.sumstats.tsv'),header=T,skip=2,sep='\t')
    tmp_stat$fam <- rep(basename(fam))
    sum_stat <- rbind(sum_stat, tmp_stat)
    mother_hom <- NULL
  }
  rm_snps<- sum_stat %>% 
    group_by(fam,Locus.ID,Col) %>%  # Each group contains one SNP (both male and female pop)
    summarise(hom_mot = sum(Pop.ID=='F' & Q.Nuc=='-'))  # If all females are homozygous -> mother is homozygous
  
  rm_snps <- rm_snps[rm_snps$hom_mot==0, 1:3]
  sum_stat <- merge(sum_stat,rm_snps,by=c("Locus.ID","Col","fam"))
}else{
  sum_stat <- read.csv(paste0(stat_path,"/batch_0.sumstats.tsv"), header=T, skip=2, sep='\t')
  sum_stat$fam <- 'all'
}


# Computing total homozygousity and number of individuals per SNP
male_stat <- sum_stat[sum_stat$Pop.ID=="M",c('Chr','BP','Obs.Hom', 'fam', 'N',"Locus.ID")]
female_stat <- sum_stat[sum_stat$Pop.ID=="F",c('Chr','BP','Obs.Hom', 'fam', 'N')]
chrom_stat <- merge(male_stat,female_stat,by=c('Chr','BP','fam'))
attach(chrom_stat)
chrom_stat$hom <- (Obs.Hom.x * N.x + Obs.Hom.y * N.y)/(N.x + N.y)
chrom_stat <- rename(chrom_stat, Nm=N.x,Nf=N.y, Male.Hom=Obs.Hom.x, Fem.Hom=Obs.Hom.y)
chrom_stat$weight <- chrom_stat$Nf + chrom_stat$Nm

# Excluding unordered contigs
chrom_stat <- chrom_stat[grep("chr.*", chrom_stat$Chr),]
chrom_stat$Chr <- droplevels(chrom_stat$Chr)


# Span of local regression has a strong effect on fit.
centrolist <- list()
centrolist[['var']] <- data.frame(Chr = levels(chrom_stat$Chr), pos = rep(0,6))
sp_range <- seq(0.1,1,0.01)
virilist <- viridis(n=length(sp_range))
colindex <- 1
par(mfrow=c(1,6))
for(chrom in levels(chrom_stat$Chr)){
  plot(x=c(),y=c(),xlim=c(0,max(chrom_stat$BP[chrom_stat$Chr==chrom])), ylim=c(-0.05,1),
       xlab=chrom,ylab="Homozygosity")
  abline(h=-0.01)
  for(sp in sp_range){
    mod <- loess(data=chrom_stat[chrom_stat$Chr==chrom,], degree=2,
                     formula=hom~BP, weights=weight, span=sp, model=T)
    points(mod$x[order(mod$x)],mod$fitted[order(mod$x)], type='l',col=alpha(virilist[colindex],0.4),lwd=1.3)
    points(mod$x[mod$fitted==min(mod$fitted)],rep(-0.01,length(mod$x[mod$fitted==min(mod$fitted)])),
           col=alpha(virilist[colindex],0.4),pch=16, cex=1.5)
    centrolist$var$pos[centrolist$var$Chr==chrom] <- mod$x[mod$fitted==min(mod$fitted)]
    colindex <- colindex+1
  }
  colindex <- 1
}

# Need to use cross-validation to select best span objectively.


loessGCV <- function (x) {
  ## Modified from code by Michael Friendly
  ## http://tolstoy.newcastle.edu.au/R/help/05/11/15899.html
  if (!(inherits(x,"loess"))) stop("Error: argument must be a loess object")
  ## extract values from loess object
  span <- x$pars$span
  n <- x$n
  traceL <- x$trace.hat
  sigma2 <- sum(resid(x)^2) / (n-1)
  gcv  <- n*sigma2 / (n-traceL)^2
  result <- list(span=span, gcv=gcv)
  result
}

estLoess <- function(model, spans = c(.1, .2)) {
  f <- function(span) {
    mod <- update(model, span = span)
    loessGCV(mod)[["gcv"]]
  }
  result <- optimize(f, spans)
  result
}

sp_range <- seq(0.2,0.5,0.05)
par(mfrow=c(3,2))
for(chrom in levels(chrom_stat$Chr)){
  plot(x=c(),y=c(),xlim=c(0,max(chrom_stat$BP[chrom_stat$Chr==chrom])), ylim=c(-0.05,1))
  abline(h=-0.01)
  mod <- loess(data = chrom_stat[chrom_stat$Chr==chrom,], hom ~ BP,
               weights = weight,degree=0)
  mod.best <- estLoess(mod)
  mod.cv <- loess(data = chrom_stat[chrom_stat$Chr==chrom,], formula=hom~BP, 
                  degree=2, weights=weight, span=mod.best$minimum, model=T)
  points(mod.cv$x[order(mod.cv$x)],mod.cv$fitted[order(mod.cv$x)], type='l')
  points(mod.cv$x[mod.cv$fitted==min(mod.cv$fitted)],rep(-0.01,length(mod.cv$x[mod.cv$fitted==min(mod.cv$fitted)])),
         pch=16, cex=1.5)
}
fix <- read.csv("../../data/assoc_mapping/prop_hom_fixed_sites.tsv", sep='\t',header=T)
fix <- fix[fix$N.Samples>0,]
fix <- fix[grep("chr.*",fix$Chr),]
fix$Chr <- droplevels(fix$Chr)
library(zoo)
HOM = zoo(fix$Prop.Hom[fix$Chr=='chr1'])
x = rollapply(HOM, width=50, by=1, FUN=mean, align='left')
plot(1:length(x),x, type='l')

HOM = zoo(chrom_stat$hom[chrom_stat$Chr=='chr1'])
x = rollapply(HOM, width=50, by=1, FUN=mean, align='left')
plot(1:length(x),x, type='l')


sp_range <- 0.75
virilist <- viridis(n=length(sp_range))
colindex <- 1
par(mfrow=c(1,6))
centrolist[['fix']] <- rep(0,6)
centrolist[['fix']] <- data.frame(pos = rep(0,6), Chr=levels(fix$Chr))
for(chrom in levels(fix$Chr)){
  plot(x=c(),y=c(),xlim=c(0,max(fix$BP[fix$Chr==chrom])), ylim=c(-0.05,1),
       xlab=chrom,ylab="Homozygosity")
  abline(h=-0.01)
  for(sp in sp_range){
    mod <- loess(data=fix[fix$Chr==chrom,], degree=2,
                 formula=Prop.Hom~BP, weights=N.Samples, span=sp, model=T)
    points(mod$x[order(mod$x)],mod$fitted[order(mod$x)], type='l',col=alpha(virilist[colindex],0.4),lwd=1.3)
    points(mod$x[mod$fitted==min(mod$fitted)],rep(-0.01,length(mod$x[mod$fitted==min(mod$fitted)])),
           col=alpha(virilist[colindex],0.4),pch=16, cex=1.5)
    centrolist$fix$pos[centrolist$fix$Chr==chrom] <- mod$x[mod$fitted==min(mod$fitted)]
    colindex <- colindex+1
  }
  colindex <- 1
}

# Plotting local regression of data with default parameters
# span = 0.75
varplot <- ggplot(chrom_stat, aes(x=BP, y=hom, weight=weight)) + facet_grid(~Chr, scales = 'free_x') + 
  geom_point(col='grey70') + geom_smooth(fill='steelblue', method='loess') + 
  geom_vline(data=centrolist$var, aes(xintercept=pos), col='red', lty=2) + theme_classic()

fixplot <- ggplot(fix, aes(x=BP, y=Prop.Hom, weight=N.Samples)) + facet_grid(~Chr, scales = 'free_x') + 
  geom_point(col='grey70') + geom_smooth(fill='steelblue', method='loess') + 
  geom_vline(data=centrolist$fix, aes(xintercept=pos), col='red', lty=2) + theme_minimal()
