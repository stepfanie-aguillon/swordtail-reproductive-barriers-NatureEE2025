# This script contains the R code to process raw data from EthoVision for CHPL samples
#
# Authors: Aguillon SM, et al.
# Year: 2024
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: TBD
# bioRxiv DOI: https://doi.org/10.1101/2024.04.16.589374
#
# Edited date: 18 Dec 2023
#
# Please cite the paper if you use these scripts
#


# load packages
library(tidyverse)
library(readxl)
library(gridExtra)

# EthoVision column names
cols <- c("trial_time", "recording_time", "x", "y", "area", "area_change", "elongation", "distance", "velocity", "zone1.1", "zone1.2", "zone2.1", "zone2.2", "result1")
##### NOTE: some of the other trials will have multiple "zone" names, so will have more columns
##### NOTE: will need to add an extra step later to combine them into a single column!


# make tibble to add data to
summary_tibble <- tibble(total_time=numeric(), trial_time=numeric(), 
                         z1_total=numeric(), z2_total=numeric(), 
                         z1_latency=numeric(), z2_latency=numeric(), 
                         z1_trial=numeric(), z2_trial=numeric(),
                         distance=numeric(), speed=numeric())

# list out all data file names to process
# these are the raw data output files from EthoVision
file_names <- Sys.glob("./data/CHPL-pheromone-trials/*.xlsx")

# frequency of tracking data collected in EthoVision
scaling_factor <- 0.2


# loop to process each file and append data
for (filename in file_names) {
  
  #reads in data file
  trial_data <- read_excel(filename, skip=36, col_names=cols, na="NA") %>%
    mutate(zone1 = zone1.1 + zone1.2, zone2 = zone2.1 + zone2.2) %>%
    ##### NOTE: can add a step here to combine columns into `zone1` and `zone2` only
    select(trial_time, x, y, zone1, zone2) %>%
    filter(trial_time <= 600)
  
  #calculate distance
  trial_data$distance <- NA
  for (i in 1:nrow(trial_data)) {
    trial_data$distance[i] <- ifelse(i==1, 0, sqrt(abs(trial_data$x[i] - trial_data$x[i-1])^2 + abs(trial_data$y[i] - trial_data$y[i-1])^2))
  }
  
  #calculate speed
  trial_data <- trial_data %>%
    mutate(speed = distance/scaling_factor)
  
  
  
  #calculate summary stats for the file
  
  #total time spent in each zone
  #validated with EthoVision summary output
  z1_total_time <- sum(trial_data$zone1, na.rm=TRUE)*scaling_factor
  z2_total_time <- sum(trial_data$zone2, na.rm=TRUE)*scaling_factor
  
  #latency to enter each zone
  #prints NA if fish never entered the zone
  #validated with EthoVision summary output
  z1_latency <- ifelse(dim(filter(trial_data, zone1==1))[1] == 0, 
                       NA, 
                       (min(select(filter(trial_data, zone1==1), trial_time))))
  z2_latency <- ifelse(dim(filter(trial_data, zone2==1))[1] == 0, 
                       NA, 
                       (min(select(filter(trial_data, zone2==1), trial_time))))
  max_latency <- max(z1_latency, z2_latency)
  
  #cut down trial df to include times AFTER the latest latency time only
  trim_trial_data <- trial_data %>%
    filter(trial_time >= max_latency)
  #calculate total time in zones AFTER the latest latency time
  z1_post_latency <- sum(trim_trial_data$zone1, na.rm=TRUE)*scaling_factor
  z2_post_latency <- sum(trim_trial_data$zone2, na.rm=TRUE)*scaling_factor
  
  #total time of the trial (max=600, but could be slightly less depending on video editing)
  total_length <- max(trial_data$trial_time)
  #calculate length of trial after the latest latency time
  trial_total = total_length - max_latency
  
  #calculate total distance traveled (cm) AFTER the latest latency time
  total_distance <- sum(trim_trial_data$distance, na.rm=TRUE)
  #calculate average speed (cm/s) AFTER the latest latency time
  avg_speed <- mean(trim_trial_data$speed, na.rm=TRUE)
  
  
  #append the  summary stats  from the individual file
  summary_tibble <- add_row(summary_tibble, total_time=total_length, trial_time=trial_total, z1_total=z1_total_time, z2_total=z2_total_time, 
                            z1_latency=z1_latency, z2_latency=z2_latency, 
                            z1_trial=z1_post_latency, z2_trial=z2_post_latency,
                            distance=total_distance, speed=avg_speed)
}



# load trial keyfile
trial_info <- read_excel("./data/CHPL-pheromone-trial-keyfile.xlsx") %>%
  arrange(file_name)

# combine the summary tibble and the keyfile
trial_summary <- cbind(trial_info, summary_tibble)





# beginning to combine the two trials for each individual

# first need to remove individuals that were unresponsive in both trials
# and any individuals that are marked to "exclude" for other reasons (e.g., issues during trials)
# then test for side bias and remove any fish that spend >80% of their time on one side

# mark unresponsive individuals
# total trial time adds up to 0 when NAs are removed (unresponsive fish have NAs, so total==0)
filter_df <- trial_summary %>%
  group_by(fish_ID) %>%
  summarize(responsiveness = ifelse(sum(trial_time, na.rm=TRUE)==0, "unresponsive", "responsive"))
trial_summary <- left_join(trial_summary, filter_df, by="fish_ID")

# remove individuals that are unresponsive or marked to exclude
trial_summary_keep <- trial_summary %>%
  filter(responsiveness=="responsive") %>%
  filter(exclude=="N")

# mark individuals that exhibit side bias and exclude
filter_df2 <- trial_summary_keep %>%
  group_by(fish_ID) %>%
  summarize(z1_bias = sum(z1_total)/sum(total_time),
            z2_bias = sum(z2_total/sum(total_time)),
            side_biased = ifelse(max(z1_bias, z2_bias)>0.80, "biased", "unbiased"))
trial_summary_keep <- left_join(trial_summary_keep, filter_df2, by="fish_ID")
trial_summary_keep <- trial_summary_keep %>%
  filter(side_biased=="unbiased")



# combining the different trials to calculate association times

# combine the preferences from the two trials for A trial lanes
#A trial lanes for trial #1
combined_trials_A1 <- trial_summary_keep %>%
  group_by(fish_ID) %>%
  filter(trial_lane=="A") %>%
  filter(trial_num==1) %>%
  summarize(cortezi1=z2_trial, birchmanni1=z1_trial, trial_time1=trial_time, 
            trial_lane=trial_lane, species=species, date=date, time=time, 
            distance1=distance, speed1=speed) %>%
  replace(is.na(.),0)
#A trial lanes for trial #2
combined_trials_A2 <- trial_summary_keep %>%
  group_by(fish_ID) %>%
  filter(trial_lane=="A") %>%
  filter(trial_num==2) %>%
  summarize(cortezi2=z1_trial, birchmanni2=z2_trial, trial_time2=trial_time,
            distance2=distance, speed2=speed) %>%
  replace(is.na(.),0)
#combine the two dfs
combined_trials_A <- left_join(combined_trials_A1, combined_trials_A2, by="fish_ID") %>%
  replace(is.na(.),0)  %>% ### NOTE: NEEDED TO ADD THIS SINCE ONE TRIAL WAS NOT REPEATED ON BOTH SIDES
  mutate(time_cortezi=cortezi1+cortezi2,
         time_birchmanni=birchmanni1+birchmanni2,
         total_time=trial_time1+trial_time2,
         prop_cortezi=time_cortezi/total_time,
         prop_birchmanni=time_birchmanni/total_time,
         prop_either=(time_cortezi+time_birchmanni)/total_time,
         SOP=(time_cortezi-time_birchmanni)/(time_cortezi+time_birchmanni),
         avg_distance=ifelse(distance1>0 & distance2>0, (distance1+distance2)/2, distance1+distance2), #addition works because one will be 0
         avg_speed=ifelse(speed1>0 & speed2>0, (speed1+speed2)/2, speed1+speed2)) %>% 
  select(date, time, fish_ID, trial_lane, species, prop_cortezi, prop_birchmanni, 
         prop_either, total_time, SOP, avg_distance, avg_speed)

# combine the preferences from the two trials for B trial lanes
#B trial lanes for trial #1
combined_trials_B1 <- trial_summary_keep %>%
  group_by(fish_ID) %>%
  filter(trial_lane=="B") %>%
  filter(trial_num==1) %>%
  summarize(cortezi1=z1_trial, birchmanni1=z2_trial, trial_time1=trial_time,
            distance1=distance, speed1=speed) %>%
  replace(is.na(.),0)
#B trial lanes for trial #2
combined_trials_B2 <- trial_summary_keep %>%
  group_by(fish_ID) %>%
  filter(trial_lane=="B") %>%
  filter(trial_num==2) %>%
  summarize(cortezi2=z2_trial, birchmanni2=z1_trial, trial_time2=trial_time, 
            trial_lane=trial_lane, species=species, date=date, time=time,
            distance2=distance, speed2=speed) %>%
  replace(is.na(.),0)
#combine the two dfs
combined_trials_B <- left_join(combined_trials_B2, combined_trials_B1, by="fish_ID") %>%
  replace(is.na(.),0)  %>% ### NOTE: NEEDED TO ADD THIS SINCE ONE TRIAL WAS NOT REPEATED ON BOTH SIDES
  mutate(time_cortezi=cortezi1+cortezi2,
         time_birchmanni=birchmanni1+birchmanni2,
         total_time=trial_time1+trial_time2,
         prop_cortezi=time_cortezi/total_time,
         prop_birchmanni=time_birchmanni/total_time,
         prop_either=(time_cortezi+time_birchmanni)/total_time,
         SOP=(time_cortezi-time_birchmanni)/(time_cortezi+time_birchmanni),
         avg_distance=ifelse(distance1>0 & distance2>0, (distance1+distance2)/2, distance1+distance2), #addition works because one will be 0
         avg_speed=ifelse(speed1>0 & speed2>0, (speed1+speed2)/2, speed1+speed2)) %>% 
  select(date, time, fish_ID, trial_lane, species, prop_cortezi, prop_birchmanni, 
         prop_either, total_time, SOP, avg_distance, avg_speed)

# combine the results from the A and B trial lanes into a single df
all_combined_trials <- rbind(combined_trials_A, combined_trials_B)







### output relevant dfs for subsequent use
#full dataset with latency, total recorded time in zones, total preference time, etc.
#each individual is separated across 2 rows
#no individuals are excluded
write.table(trial_summary, "./CHPL-pheromone-overall-trial-summary.txt", quote=FALSE, sep="\t", row.names=FALSE)

#filtered dataset with latency, total recorded time in zones, total preference time, etc.
#each individual is separated across 2 rows
write.table(trial_summary_keep, "./CHPL-pheromone-filtered-trial-summary.txt", quote=FALSE, sep="\t", row.names=FALSE)

#filtered dataset that includes overall cue "preference" for each individual
#excludes individuals that don't pass conditionals, etc.
write.table(all_combined_trials, "./CHPL-pheromone-trial-preference.txt", quote=FALSE, sep="\t", row.names=FALSE)


