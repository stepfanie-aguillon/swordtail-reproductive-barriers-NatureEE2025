# This script contains the R code to analyze video observations
#
# Authors: Aguillon SM, et al.
# Year: 2024
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: TBD
# bioRxiv DOI: https://doi.org/10.1101/2024.04.16.589374
#
# Edited date: 5 Dec 2024
#
# Please cite the paper if you use these scripts
#

# load packages
library(tidyverse)

# load data
obs <- read_csv("./data/video-observations-summary.csv", col_names=TRUE) %>%
  filter(Duration_min > 1) %>%
  select(Video:`M together`) %>%
  mutate(total_time_sec = Duration_min*60 + Duration_sec) # calculate total video time in seconds

# total video observation time
sum(obs$total_time_sec/60)  #247.9 minutes (4.131667 hours)


# total number of sworded males observed
sum(obs$Sworded) #70

# total number of unsworded males observed
sum(obs$Unsworded) #386

# total number of females 
sum(obs$Females) #756


