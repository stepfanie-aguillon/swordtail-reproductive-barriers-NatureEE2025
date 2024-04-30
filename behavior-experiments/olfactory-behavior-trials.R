# This script contains the R code to analyze EthoVision behavior results
#
# Authors: Aguillon SM, et al.
# Year: 2024
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: TBD
# bioRxiv DOI: https://doi.org/10.1101/2024.04.16.589374
#
# Edited date: 22 Jan 2024
#
# Please cite the paper if you use these scripts

# load packages
library(tidyverse)
library(cowplot)
library(ggrepel)

# plot colors
corcol <- "#359D6F"
bircol <- "#3A6BC6"


### load data ###

# olfactory trials using allopatric females
allopatric_summary <- read_tsv("./pheromone-filtered-trial-summary.txt", col_names=TRUE)
allopatric_preference <- read_tsv("./pheromone-trial-preference.txt", col_names=TRUE) %>%
  mutate(prop_neutral = 1 - prop_either)

# olfactory trials using females from Chapulhuacanito
CHPL_summary <- read_tsv("./CHPL-pheromone-filtered-trial-summary.txt", col_names=TRUE)
CHPL_preference <- read_tsv("./CHPL-pheromone-trial-preference.txt", col_names=TRUE) 
# females of unknown ancestry, and later sequenced
CHPL_ancestry <- read_tsv("./data/CHPL-contemporary-ancestry_behaviortrial-subset.txt", col_names=TRUE)

CHPL_preference_2 <- left_join(CHPL_preference, CHPL_ancestry, by="fish_ID") %>%
  mutate(ancestry_cluster = ifelse(hybrid_index < 0.10, "birchmanni", "cortezi-like"),
         prop_neutral = 1 - prop_either)



##### behavior trial results: allopatric samples, SOP #####
allopatric_SOP <- ggplot(data=allopatric_preference, aes(x=species, y=SOP)) +
  geom_hline(yintercept=0,color="dark gray") + 
  geom_boxplot(aes(fill=species), outlier.shape=NA) +
  geom_jitter(color="gray", width=0.1) +
  scale_fill_manual(values=c(bircol, corcol)) + 
  theme_bw() +
  xlab("Species") +
  ylab("Relative preference for cortezi") +
  ylim(c(-1.0,1.1)) + 
  scale_x_discrete(labels=c("birchmanni", "cortezi")) +
  theme(legend.position="none", axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"),strip.text.x=element_text(size=12,face="bold"))

# statistically test results for SOP (Strength of Preference)
## allopatric trials
wilcox.test(filter(allopatric_preference, species=="birchmanni")$SOP, mu=0, alternative="less")
#Wilcoxon signed rank test with continuity correction
#
#data:  filter(allopatric_preference, species == "birchmanni")$SOP
#V = 83, p-value = 0.7886
#alternative hypothesis: true location is less than 0
wilcox.test(filter(allopatric_preference, species=="cortezi")$SOP, mu=0, alternative="great")
#data:  filter(allopatric_preference, species == "cortezi")$SOP
#V = 94, p-value = 0.2105
#alternative hypothesis: true location is greater than 0


##### behavior trial results: allopatric samples, neutral zone #####
allopatric_neutral <- ggplot(data=allopatric_preference, aes(x=species, y=prop_neutral)) +
  geom_hline(yintercept=1/3,color="dark gray") + 
  geom_boxplot(aes(fill=species), outlier.shape=NA) +
  geom_jitter(color="gray", width=0.1) +
  scale_fill_manual(values=c(bircol, corcol)) + 
  theme_bw() +
  xlab("Species") +
  ylab("Proportion time in neutral zone") +
  ylim(c(-0.05,1.05)) + 
  scale_x_discrete(labels=c("birchmanni", "cortezi")) +
  theme(legend.position="none", axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"),strip.text.x=element_text(size=12,face="bold"))

# statistically test results for neutral zone
## allopatric trials
wilcox.test(filter(allopatric_preference, species=="birchmanni")$prop_neutral, mu=1/3)
#Wilcoxon signed rank exact test
#
#data:  filter(allopatric_preference, species == "birchmanni")$prop_neutral
#V = 99, p-value = 0.1167
#alternative hypothesis: true location is not equal to 0.3333333
wilcox.test(filter(allopatric_preference, species=="cortezi")$prop_neutral, mu=1/3)
#Wilcoxon signed rank exact test
#
#data:  filter(allopatric_preference, species == "cortezi")$prop_neutral
#V = 14, p-value = 0.001678
#alternative hypothesis: true location is not equal to 0.3333333

wilcox.test(filter(allopatric_preference, species=="birchmanni")$prop_neutral, filter(allopatric_preference, species=="cortezi")$prop_neutral)
#Wilcoxon rank sum exact test
#
#data:  filter(allopatric_preference, species == "birchmanni")$prop_neutral and filter(allopatric_preference, species == "cortezi")$prop_neutral
#W = 221, p-value = 0.001634
#alternative hypothesis: true location shift is not equal to 0






##### behavior trial results: CHPL samples #####
CHPL_SOP <- ggplot(data=CHPL_preference_2, aes(x=ancestry_cluster, y=SOP)) +
  geom_hline(yintercept=0,color="dark gray") + 
  geom_boxplot(aes(fill=ancestry_cluster), outlier.shape=NA) +
  geom_jitter(color="gray", width=0.1) +
  scale_fill_manual(values=c(bircol, corcol)) + 
  theme_bw() +
  xlab("Ancestry cluster") +
  ylab("Relative preference for cortezi") +
  ylim(c(-1.0,1.1)) + 
  scale_x_discrete(labels=c("birchmanni", "cortezi-like")) +
  theme(axis.title.y=element_blank(), legend.position="none", axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"),strip.text.x=element_text(size=12,face="bold"))

# statistically test results for SOP (Strength of Preference)
## allopatric trials
wilcox.test(filter(CHPL_preference_2, ancestry_cluster=="birchmanni")$SOP, mu=0, alternative="less")
#Wilcoxon signed rank exact test
#
#data:  filter(CHPL_preference_2, ancestry_cluster == "birchmanni")$SOP
#V = 60, p-value = 0.511
#alternative hypothesis: true location is less than 0
wilcox.test(filter(CHPL_preference_2, ancestry_cluster=="cortezi-like")$SOP, mu=0, alternative="great")
#Wilcoxon signed rank exact test
#
#data:  filter(CHPL_preference_2, ancestry_cluster == "cortezi-like")$SOP
#V = 49, p-value = 0.5961
#alternative hypothesis: true location is greater than 0


##### behavior trial results: CHPL samples, neutral zone #####
CHPL_neutral <- ggplot(data=CHPL_preference_2, aes(x=ancestry_cluster, y=prop_neutral)) +
  geom_hline(yintercept=1/3,color="dark gray") + 
  geom_boxplot(aes(fill=ancestry_cluster), outlier.shape=NA) +
  geom_jitter(color="gray", width=0.1) +
  scale_fill_manual(values=c(bircol, corcol)) + 
  theme_bw() +
  xlab("Ancestry cluster") +
  ylab("Proportion time in neutral zone") +
  ylim(c(-0.05,1.05)) + 
  scale_x_discrete(labels=c("birchmanni", "cortezi-like")) +
  theme(axis.title.y=element_blank(), legend.position="none", axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"),strip.text.x=element_text(size=12,face="bold"))

# statistically test results for neutral zone
## allopatric trials
wilcox.test(filter(CHPL_preference_2, ancestry_cluster=="birchmanni")$prop_neutral, mu=1/3)
#Wilcoxon signed rank exact test
#
#data:  filter(CHPL_preference_2, ancestry_cluster == "birchmanni")$prop_neutral
#V = 79, p-value = 0.3028
#alternative hypothesis: true location is not equal to 0.3333333
wilcox.test(filter(CHPL_preference_2, ancestry_cluster=="cortezi-like")$prop_neutral, mu=1/3)
#Wilcoxon signed rank exact test
#
#data:  filter(CHPL_preference_2, ancestry_cluster == "cortezi-like")$prop_neutral
#V = 0, p-value = 0.0001221
#alternative hypothesis: true location is not equal to 0.3333333

wilcox.test(filter(CHPL_preference_2, ancestry_cluster=="birchmanni")$prop_neutral, filter(CHPL_preference_2, ancestry_cluster=="cortezi-like")$prop_neutral)
#Wilcoxon rank sum exact test
#
#data:  filter(CHPL_preference_2, ancestry_cluster == "birchmanni")$prop_neutral and filter(CHPL_preference_2, ancestry_cluster == "cortezi-like")$prop_neutral
#W = 197, p-value = 9.619e-06
#alternative hypothesis: true location shift is not equal to 0





#### plot the combined results for the supplement, SOP ####
plot_grid(allopatric_SOP, CHPL_SOP, align="h", nrow=1)
# export as 6.5" x 4" PDF
# in Illustrator: italicized species names, fixed plot titles, added panel letters


#### plot the combined results for the supplement, neutral zone ####
plot_grid(allopatric_neutral, CHPL_neutral, align="h", nrow=1)
# export as 6.5" x 4" PDF
# in Illustrator: italicized species names, fixed plot titles, added panel letters


