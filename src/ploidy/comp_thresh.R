# This script produces visualisations to compare the different thresholds of ploidy separation.
# It takes the summary tables output by the different thresholds used in haplo_males.py.
# Cyril Matthey-Doret
# 26.05.2017

suppressMessages(library(tidyverse))

#in_folder <- '../../data/ploidy/thresholds/'
in_folder <- commandArgs(trailingOnly = TRUE)  # Folder containing input tables

#' ploidy_plot
#' Produces one barplot per families representing the inbreeding coefficient
#' of all individuals, along with their state (mother, daughter, haploid son,
#' diploid son)
#' @param ploid_tbl input table containing all individuals 
#' @param thresh filename of the table, showing threshold multiplier and transformation
#'
#' @return Nothing, saves a plot into a folder
ploidy_plot <- function(ploid_tbl,thresh, type='bar'){
  
  group_stats <- ploid_tbl %>%  # Computing summary stats for each family
    filter(Generation == 'F4' & Sex == 'F') %>%  # Computing from daughters only
    group_by(Family) %>%
    summarise_at(.cols=c("Fis"),.funs = c("mean","sd")) %>%
    mutate(state = "Daughters") %>%  # These are all daughters (used for plotting)
    mutate(mid_x = n()/2)  # Graphical parameter for plots
  group_stats$sd[is.na(group_stats$sd)] <- 0
  # Producing barplot faceted by family, with mean and standard deviation displayed
  # as segments. Colors represent individuals' states.
  if(type=='bar'){
    gg <- ggplot(data = ploid_tbl, aes(x = factor(Fis), y = Fis,fill = state))+
      geom_bar(stat = 'identity') +
      geom_hline(data = group_stats, aes(yintercept = mean)) + 
      geom_segment(data = group_stats,aes(x=mid_x,xend=mid_x,y=mean-sd,yend=mean+sd)) +
      theme(axis.text.x = element_blank()) + facet_wrap(~Family,drop=T, scale='free') +
      xlab("Individuals") + ylab("Inbreeding coefficient") + ggtitle(thresh)
    out_folder <- 'barplots'
  } else{
    gg <- ggplot(data = ploid_tbl, aes(x = Fis))+
      geom_density(data = ploid_tbl[ploid_tbl$state!="Mothers",], alpha = 0.3,aes(fill=Sex)) +  # Daughters vs Sons as density area
      #geom_density(alpha = 0.1, col="green") +  # All individuals in family showed by green density curve
      geom_point(data = ploid_tbl[ploid_tbl$state=="Mothers",], aes(y=0, col=state)) +  # Showing mother as a dot
      geom_vline(data = group_stats, aes(xintercept = (mean + 2* sd)), col='red', lty=2) +
      theme(axis.text.x = element_blank()) + facet_wrap(~Family,drop=T, scale='free') +
      xlab("Inbreeding coefficient") + ylab("Density") + ggtitle(paste("Threshold ", thresh, sep=" "))
    out_folder <- 'density'
  }
  pdf(paste0('data/ploidy/plots/', out_folder,'/',thresh,'.pdf'))  # Opening pdf connection
  print(gg)  # Saving plot
  dev.off()  # Closing connection
}

for(in_file in list.files(in_folder)){  # Iterating over input lists
  in_data <- read_tsv(paste0(in_folder,in_file),col_names = T)  # Reading files sequentially
  
  # Preparing data:
  in_data %<>%   # Erasing original variable
    mutate(state=paste0(Generation,Sex,Ploidy)) %>%  # New column containing state information
    rename(Fis=F)  # Renaming inbreeding coefficient column (F is colliding with 'FALSE')
  # Changing individual states by Human-friendly names
  in_data$state <-recode_factor(in_data$state, "F3FD" = "Mothers",
                         "F4FD" = "Daughters",
                         "F4MH" = "1N Sons", 
                         "F4MD" = "2N Sons")
  ploidy_plot(in_data,basename(in_file))  # Calling plotting function
  ploidy_plot(in_data,basename(in_file), type='density')
}
