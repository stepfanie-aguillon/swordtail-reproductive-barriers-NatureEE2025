# This script contains the R code to assess ancestryinfer priors
#
# Authors: Aguillon SM, et al.
# Year: 2025
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: Nature Ecology and Evolution
# DOI: https://doi.org/10.1038/s41559-025-02669-9
#
# Edited date: 14 Dec 2023
#
# Please cite the paper if you use these scripts
#


# load packages
library(tidyverse)
library(cowplot)

# load data and process for summary
xbir <- read_tsv("./data/testing-priors-xbir.txt", col_names=TRUE) %>%
  mutate(diff=abs(hybrid_index_N - hybrid_index_U)*100, cluster="xbir")
xcor <- read_tsv("./data/testing-priors-xcor.txt", col_names=TRUE) %>%
  mutate(diff=abs(hybrid_index_N - hybrid_index_U)*100, cluster="xcor")

# calculate absolute value of difference between uniform priors and cluster-specific priors
# calculate mean and standard deviation of these differences
combined_df <- bind_rows(xbir, xcor) 

combined_df %>%
  group_by(cluster) %>%
  summarize(n = n(), 
            mean_diff_percent = mean(diff),
            sd_diff_percent = sd(diff))
## A tibble: 2 × 4
#cluster     n mean_diff_percent sd_diff_percent
#<chr>   <int>             <dbl>           <dbl>
#1 xbir      102            0.0900          0.0244
#2 xcor       69            0.495           0.514 


# calculate overall mean and standard deviation for all samples
mean(combined_df$diff) #0.2535958
sd(combined_df$diff) #0.3817439

# test for differences between the two clusters
wilcox.test(xcor$diff, xbir$diff)
#Wilcoxon rank sum test with continuity correction
#
#data:  xcor$diff and xbir$diff
#W = 6614, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0
wilcox.test(xcor$diff, xbir$diff)$p.value
# 1.976515e-22

# test for differences in each cluster from 0

#xcor-admixed cluster
wilcox.test(xcor$diff, mu=0)
#Wilcoxon signed rank test with continuity correction
#
#data:  xcor$diff
#V = 2415, p-value = 5.331e-13
#alternative hypothesis: true location is not equal to 0

#xbir-admixed cluster
wilcox.test(xbir$diff, mu=0)
#Wilcoxon signed rank test with continuity correction
#
#data:  xbir$diff
#V = 5253, p-value < 2.2e-16
#alternative hypothesis: true location is not equal to 0
wilcox.test(xbir$diff, mu=0)$p.value
#1.850183e-18 


#plot for xbir cluster
xbir_plot <- ggplot(data=xbir) + 
  geom_point(aes(x=hybrid_index_U, y=hybrid_index_N), color="#3A6BC6") +
  geom_abline(intercept=0) + 
  ylab("Hybrid index (admixture prior = 0.99)") +
  theme_bw() + 
  theme(axis.title.x=element_blank(), axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"), strip.text.x=element_text(size=12,face="bold"))


#plot for admixed-xcor cluster
xcor_plot <- ggplot(data=xcor) + 
  geom_point(aes(x=hybrid_index_U, y=hybrid_index_N), color="#359D6F") + 
  geom_abline(intercept=0) +
  xlab("Hybrid index (admixture prior = 0.50)") +
  ylab("Hybrid index (admixture prior = 0.75)") +
  theme_bw() + 
  theme(axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"), strip.text.x=element_text(size=12,face="bold"))


#### plot the combined results for the supplement ####
plot_grid(xbir_plot, xcor_plot, align="v", nrow=2)
# export as 5" x 8" PDF
# in Illustrator: added panel labels and cluster information 
