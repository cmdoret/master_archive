# This script produces visualisations to compare the different thresholds of ploidy separation.
# It takes the summary tables output by the different thresholds used in haplo_males.py.
# Cyril Matthey-Doret
# 26.05.2017

library(tidyverse)

#in_folder <- '../../data/ploidy/'
in_folder <- commandArgs(trailingOnly = TRUE)  # Folder containing input tables

#' ploidy_plot
#' Produces one barplot per families representing the inbreeding coefficient
#' of all individuals, along with their state (mother, daughter, haploid son,
#' diploid son)
#' @param ploid_tbl input table containing all individuals 
#' @param thresh filename of the table, showing threshold multiplier and transformation
#'
#' @return Nothing, saves a plot into a folder
ploidy_plot <- function(ploid_tbl,thresh){
  
  group_stats <- ploid_tbl %>%  # Computing summary stats for each family
    filter(Generation == 'F4' & Sex == 'F') %>%  # Computing from daughters only
    group_by(Family) %>%
    summarise_at(.cols=c("Fis"),.funs = c("mean","sd")) %>%
    mutate(state = "Daughters") %>%  # These are all daughters (used for plotting)
    mutate(mid_x = n()/2)  # Graphical parameter for plots
  
  # Producing barplot faceted by family, with mean and standard deviation displayed
  # as segments. Colors represent individuals' states.
  gg <- ggplot(data = ploid_tbl, aes(x = factor(Fis), y = Fis,fill = state))+
    geom_bar(stat = 'identity') +
    geom_hline(data = group_stats, aes(yintercept = mean)) + 
    geom_segment(data = group_stats,aes(x=mid_x,xend=mid_x,y=mean-sd,yend=mean+sd)) +
    theme(axis.text.x = element_blank()) + facet_wrap(~Family,drop=T, scale='free') +
    xlab("Individuals") + ylab("Inbreeding coefficient") + ggtitle(thresh)
  
  pdf(paste0('data/ploidy/plots/',thresh,'.pdf'))  # Opening pdf connection
  print(gg)  # Saving plot
  dev.off()  # Closing connection
}

haplo_prop <- tibble(threshold = numeric(0),
                     Daughters = numeric(0),
                     Mothers = numeric(0)
                     )
for(in_file in list.files(in_folder)){  # Iterating over input lists
  in_data <- read_tsv(paste0(in_folder,in_file),col_names = T)  # Reading files sequentially
  
  # Preparing data:
  in_data %<>%   # Erasing original variable
    mutate(state=paste0(Generation,Sex,Ploidy)) %>%  # New column containing state information
    rename(Fis=F)  # Renaming inbreeding coefficient column (F is colliding with 'FALSE')
  # Changing individual states by Human-friendly names
  in_data$state <-recode_factor(in_data$state, "F3FD" = "Mothers",  #
                         "F4FD" = "Daughters",
                         "F4MH" = "1N Sons", 
                         "F4MD" = "2N Sons")
  ploidy_plot(in_data,basename(in_file))  # Calling plotting function
}
