# This script contains the R code to visualize sperm results
#
# Authors: Aguillon SM, et al.
# Year: 2025
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: Nature Ecology and Evolution
# DOI: https://doi.org/10.1038/s41559-025-02669-9
#
# Edited date: 31 October 2024
#
# Please cite the paper if you use these scripts
#

# load packages
library(tidyverse)
library(readxl)
library(cowplot)


#### supplementary figures for sperm morphology data ####

# load data
cols = c("ID", "species", "obj", "sperm_number", "head_l", "head_w", "mid_l")
sperm_df <- read_excel("/Users/stepfanie/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Xcortezi_Xbirchmanni_reproductive_barriers/Data/Sperm_imaging/Data/corbir_sperm_morpho.xlsx",
                       col_names=cols, na="NA", skip=1) %>%
  filter(obj == 40) %>%
  select(ID, species, sperm_number, head_l, head_w, mid_l) %>%
  mutate(head_shape = head_l/head_w)

# calculate averages across the 10 observations / individual
sperm_means <- sperm_df %>%
  group_by(ID) %>%
  reframe(species = species, avg_head_l = mean(head_l), avg_head_w = mean(head_w), avg_mid_l = mean(mid_l), avg_head_shape = mean(head_shape)) %>%
  distinct()


### pairwise statistics - morphology

# set up the four groups
xbir <- sperm_means %>%
  filter(species=="bir")
xcor <- sperm_means %>%
  filter(species=="cor")
f1_f2 <- sperm_means %>%
  filter(species=="f1" | species=="f2")
xbir_f1_f2 <- sperm_means %>%
  filter(species=="f1" | species=="f2" | species=="bir")
xcor_f1_f2 <- sperm_means %>%
  filter(species=="f1" | species=="f2" | species=="cor")

# pairwise wilcoxons

# head length
wilcox.test(xbir$avg_head_l, xcor_f1_f2$avg_head_l, alternative="two.sided")
#data:  xbir$avg_head_l and xcor_f1_f2$avg_head_l
#W = 48, p-value = 0.001099
#alternative hypothesis: true location shift is not equal to 0


# head width
wilcox.test(xbir$avg_head_w, f1_f2$avg_head_w, alternative="two.sided")
#data:  xbir$avg_head_w and f1_f2$avg_head_w
#W = 30, p-value = 0.01616
#alternative hypothesis: true location shift is not equal to 0


# midpiece length
wilcox.test(xcor$avg_mid_l, xbir_f1_f2$avg_mid_l, alternative="two.sided")
#data:  xcor$avg_mid_l and xbir_f1_f2$avg_mid_l
#W = 48, p-value = 0.001099
#alternative hypothesis: true location shift is not equal to 0



# plot colors
corcol <- "#359D6F"
bircol <- "#3A6BC6"

head_l_plot <- ggplot(data=sperm_means, aes(x=species, y=avg_head_l)) +
  geom_boxplot(aes(fill=species), outlier.shape=NA) + 
  geom_jitter(color="gray", width=0.1) +
  scale_fill_manual(values=c(bircol, corcol, "#595959", "#595959")) +
  theme_bw() + 
  ylab("Head length (um)") + 
  scale_x_discrete(labels=c("X. birchmanni" , "X. cortezi", "F1", "F2")) +
  theme(legend.position="none", axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"),strip.text.x=element_text(size=10,face="bold")) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())

head_w_plot <- ggplot(data=sperm_means, aes(x=species, y=avg_head_w)) +
  geom_boxplot(aes(fill=species), outlier.shape=NA) + 
  geom_jitter(color="gray", width=0.1) +
  scale_fill_manual(values=c(bircol, corcol, "#595959", "#595959")) +
  theme_bw() + 
  ylab("Head width (um)") + 
  scale_x_discrete(labels=c("X. birchmanni" , "X. cortezi", "F1", "F2")) +
  theme(legend.position="none", axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"),strip.text.x=element_text(size=10,face="bold")) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())

head_shape_plot <- ggplot(data=sperm_means, aes(x=species, y=avg_head_shape)) +
  geom_boxplot(aes(fill=species), outlier.shape=NA) + 
  geom_jitter(color="gray", width=0.1) +
  scale_fill_manual(values=c(bircol, corcol, "#595959", "#595959")) +
  theme_bw() + 
  ylab("Head shape (length/width)") + 
  scale_x_discrete(labels=c("X. birchmanni" , "X. cortezi", "F1", "F2")) +
  theme(legend.position="none", axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"),strip.text.x=element_text(size=10,face="bold")) +
  theme(axis.text.x=element_text(size=10), axis.title.x=element_blank())

mid_l_plot <- ggplot(data=sperm_means, aes(x=species, y=avg_mid_l)) +
  geom_boxplot(aes(fill=species), outlier.shape=NA) + 
  geom_jitter(color="gray", width=0.1) +
  scale_fill_manual(values=c(bircol, corcol, "#595959", "#595959")) +
  theme_bw() + 
  ylab("Midpiece length (um)") + 
  scale_x_discrete(labels=c("X. birchmanni" , "X. cortezi", "F1", "F2")) +
  theme(legend.position="none", axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"),strip.text.x=element_text(size=10,face="bold")) +
  theme(axis.text.x=element_text(size=10), axis.title.x=element_blank())

#combine plots
plot_grid(head_l_plot, head_w_plot, head_shape_plot, mid_l_plot, align="hv", nrow=2, labels="AUTO", label_size=12)
# export as 6" tall x 6.5" wide
# in Illustrator: italicize species names, fix spacing slightly, update `u` to `mu symbol`



#### supplementary figures for sperm motility data ####

# load data
sperm_motility_df <- read_csv("/Users/stepfanie/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Xcortezi_Xbirchmanni_reproductive_barriers/Data/Sperm_imaging/Data/mean_motilities.csv", 
                              col_names=TRUE)

### pairwise statistics - motility

# set up the four groups
xbir <- sperm_motility_df %>%
  filter(geno=="bir")
xcor <- sperm_motility_df %>%
  filter(geno=="cor")
f1_f2 <- sperm_motility_df %>%
  filter(geno=="f1" | geno=="f2")
xbir_f1_f2 <- sperm_motility_df %>%
  filter(geno=="f1" | geno=="f2" | geno=="bir")
xcor_f1_f2 <- sperm_motility_df %>%
  filter(geno=="f1" | geno=="f2" | geno=="cor")


# pairwise comparisons

wilcox.test(xbir$`VCL Mean (um/s)`, xcor_f1_f2$`VCL Mean (um/s)`, alternative="two.sided")
#data:  xbir$`VCL Mean (um/s)` and xcor_f1_f2$`VCL Mean (um/s)`
#W = 46, p-value = 0.004396
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(xbir$`VAP Mean (um/s)`, xcor_f1_f2$`VAP Mean (um/s)`, alternative="two.sided")
#data:  xbir$`VAP Mean (um/s)` and xcor_f1_f2$`VAP Mean (um/s)`
#W = 48, p-value = 0.001099
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(xbir$`STR Mean`, xcor_f1_f2$`STR Mean`, alternative="two.sided")
#data:  xbir$`STR Mean` and xcor_f1_f2$`STR Mean`
#W = 8, p-value = 0.05824
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(xbir$`VSL Mean (um/s)`, xcor_f1_f2$`VSL Mean (um/s)`, alternative="two.sided")
#data:  xbir$`VSL Mean (um/s)` and xcor_f1_f2$`VSL Mean (um/s)`
#W = 48, p-value = 0.001099
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(xbir$`Progressive Motility (%)`, xcor_f1_f2$`Progressive Motility (%)`, alternative="two.sided")
#data:  xbir$`Progressive Motility (%)` and xcor_f1_f2$`Progressive Motility (%)`
#W = 48, p-value = 0.004346
#alternative hypothesis: true location shift is not equal to 0


VAP_plot <- ggplot(data=sperm_motility_df, aes(x=geno, y=`VAP Mean (um/s)`)) +
  geom_boxplot(aes(fill=geno), outlier.shape=NA) + 
  geom_jitter(color="gray", width=0.1) +
  scale_fill_manual(values=c(bircol, corcol, "#595959", "#595959")) +
  theme_bw() + 
  ylab("Avg path velocity (um/s)") + 
  scale_x_discrete(labels=c("X. birchmanni" , "X. cortezi", "F1", "F2")) +
  theme(legend.position="none", axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"),strip.text.x=element_text(size=10,face="bold")) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())

STR_plot <- ggplot(data=sperm_motility_df, aes(x=geno, y=`STR Mean`)) +
  geom_boxplot(aes(fill=geno), outlier.shape=NA) + 
  geom_jitter(color="gray", width=0.1) +
  scale_fill_manual(values=c(bircol, corcol, "#595959", "#595959")) +
  theme_bw() + 
  ylab("Path straightness") + 
  scale_x_discrete(labels=c("X. birchmanni" , "X. cortezi", "F1", "F2")) +
  theme(legend.position="none", axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"),strip.text.x=element_text(size=10,face="bold")) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())

PROGMOT_plot <- ggplot(data=sperm_motility_df, aes(x=geno, y=`Progressive Motility (%)`)) +
  geom_boxplot(aes(fill=geno), outlier.shape=NA) + 
  geom_jitter(color="gray", width=0.1) +
  scale_fill_manual(values=c(bircol, corcol, "#595959", "#595959")) +
  theme_bw() + 
  ylab("Progressive motility (%)") + 
  scale_x_discrete(labels=c("X. birchmanni" , "X. cortezi", "F1", "F2")) +
  theme(legend.position="none", axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"),strip.text.x=element_text(size=10,face="bold")) +
  theme(axis.text.x=element_text(size=10), axis.title.x=element_blank())

#combine plots
plot_grid(VAP_plot, STR_plot, PROGMOT_plot, align="h", nrow=3, labels="AUTO", label_size=12)
# export as 7.5" tall x 4" wide
# in Illustrator: italicize species names, update `u` to `mu symbol`
