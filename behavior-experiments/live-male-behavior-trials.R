# This script contains the R code to work with behavior data from Carla Gutierrez
#
# Authors: Aguillon SM, et al.
# Year: 2025
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: Nature Ecology and Evolution
# DOI: https://doi.org/10.1038/s41559-025-02669-9
#
# Edited date: 24 Oct 2024
#
# Please cite the paper if you use these scripts
#


# load packages
library(tidyverse)
library(cowplot)
library(readxl)

# plot colors
corcol <- "#359D6F"
bircol <- "#3A6BC6"



# load datafile
df <- read_excel("./data/live-male-behavior-data.xlsx", sheet=1, col_names=TRUE, na="NA", trim_ws=TRUE)

# trial length
trial_time <- 600*2

# remove samples where males were marked as having FGS after trials
# remove samples from Pumphouse, which is now known to be a hybrid population
df_keep <- df %>%
  filter(orig_removed=="N")

# test for side bias and remove any fish that spend >80% of their time on one side
biased_df <- df_keep %>%
  group_by(new_ID) %>%
  summarize(pair1_bias1 = (Xbir_t1+Xcor_t2)/trial_time,
            pair1_bias2 = (Xcor_t1+Xbir_t2)/trial_time,
            pair2_bias1 = (Xbir_t3+Xcor_t4)/trial_time, 
            pair2_bias2 = (Xcor_t3+Xbir_t4)/trial_time,
            pair1_side_biased = ifelse(max(pair1_bias1, pair1_bias2)>0.800, "biased", "unbiased"),
            pair2_side_biased = ifelse(max(pair2_bias1, pair2_bias2)>0.800, "biased", "unbiased"))

joined_df <- left_join(df_keep, biased_df, by="new_ID")


# separate the two pairs of trials from each other, remove samples with side bias
pair1_trials <- select(joined_df, new_ID:male_pair_ID, Xbir_t1:Xcor_t2, trial_type, pair1_bias1:pair1_bias2, pair1_side_biased) %>%
  filter(pair1_side_biased=="unbiased")
pair2_trials <- select(joined_df, new_ID:male_pair_ID, Xbir_t3:Xcor_t4, trial_type, pair2_bias1:pair2_bias2, pair2_side_biased) %>%
  filter(pair2_side_biased=="unbiased")

# calculate association times
pair1_trials <- pair1_trials %>%
  mutate(Xcor_pref = Xcor_t1 + Xcor_t2 - Xbir_t1 - Xbir_t2,
         SOP = Xcor_pref/(Xcor_t1 + Xcor_t2 + Xbir_t1 + Xbir_t2),
         trial_num = "first") %>%
  select(new_ID, species, trial_type, Xcor_pref, trial_num, SOP)

pair2_trials <- pair2_trials %>%
  mutate(Xcor_pref = Xcor_t3 + Xcor_t4 - Xbir_t3 - Xbir_t4,
         SOP = Xcor_pref/(Xcor_t3 + Xcor_t4 + Xbir_t3 + Xbir_t4),
         trial_num = "second") %>%
  select(new_ID, species, trial_type, Xcor_pref, trial_num, SOP)

full_assoc <- bind_rows(pair1_trials, pair2_trials) %>%
  mutate(trial_type=ifelse(trial_type=="vis+odor","visual+olfactory","visual"))






##### behavior trial results for in-text figure (Figure 2B) #####
ggplot(data=full_assoc, aes(x=species, y=SOP)) +
  facet_grid(~factor(trial_type, levels=c("visual","visual+olfactory"))) +
  geom_hline(yintercept=0,color="dark gray") + 
  geom_boxplot(aes(fill=species), outlier.shape=NA) +
  geom_jitter(color="gray", width=0.1) +
  scale_fill_manual(values=c(bircol, corcol)) + 
  theme_bw() +
  xlab("Species") +
  ylab("Relative preference for cortezi") +
  ylim(c(-1.0,1.1)) + 
  scale_x_discrete(labels=c("birchmanni", "cortezi")) +
  theme(legend.position="none", axis.title=element_text(face="bold",size=10), axis.text=element_text(size=8,color="black"),strip.text.x=element_text(size=10,face="bold"))
# export as 3.5" x 3" PDF
# in Illustrator: italicized species names, fixed plot titles, added significance asterisks
# March 2024 note: ended up combining this with figure on mother/embryo sequencing!



#### summarize numbers of individuals included in each analysis
full_assoc <- full_assoc %>%
  group_by(trial_type, species)

summarize(full_assoc,
          n=n())
## A tibble: 4 × 3
## Groups:   trial_type [2]
#trial_type       species        n
#<chr>            <chr>      <int>
#1 visual           birchmanni    41
#2 visual           cortezi       18
#3 visual+olfactory birchmanni    37
#4 visual+olfactory cortezi       30






# statistically test results for SOP (Strength of Preference)
## visual trials
wilcox.test(filter(full_assoc, species=="birchmanni" & trial_type=="visual")$SOP, mu=0, alternative="less")
#Wilcoxon signed rank exact test
#
#data:  filter(full_assoc, species == "birchmanni" & trial_type == "visual")$SOP
#V = 416, p-value = 0.4288
#alternative hypothesis: true location is less than 0
wilcox.test(filter(full_assoc, species=="cortezi" & trial_type=="visual")$SOP, mu=0, alternative="great")
#Wilcoxon signed rank test with continuity correction
#
#data:  filter(full_assoc, species == "cortezi" & trial_type == "visual")$SOP
#V = 144, p-value = 0.005765
#alternative hypothesis: true location is greater than 0

## visual + olfactory trials
wilcox.test(filter(full_assoc, species=="birchmanni" & trial_type=="visual+olfactory")$SOP, mu=0, alternative="less")
#Wilcoxon signed rank test with continuity correction
#
#data:  filter(full_assoc, species == "birchmanni" & trial_type == "visual+olfactory")$SOP
#V = 347.5, p-value = 0.5932
#alternative hypothesis: true location is less than 0
wilcox.test(filter(full_assoc, species=="cortezi" & trial_type=="visual+olfactory")$SOP, mu=0, alternative="great")
#Wilcoxon signed rank test with continuity correction
#
#data:  filter(full_assoc, species == "cortezi" & trial_type == "visual+olfactory")$SOP
#V = 196, p-value = 0.6829
#alternative hypothesis: true location is greater than 0
