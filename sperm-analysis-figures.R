# This script contains the R code to visualize sperm results
#
# Authors: Aguillon SM, et al.
# Year: 2024
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: TBD
# bioRxiv DOI: TBD
#
# Edited date: 13 April 2024
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
sperm_df <- read_excel("/Users/aguillon/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Xcortezi_Xbirchmanni_reproductive_barriers/Data/Sperm_imaging/Data/corbir_sperm_morpho.xlsx",
                       col_names=cols, na="NA", skip=1) %>%
  filter(obj == 40) %>%
  select(ID, species, sperm_number, head_l, head_w, mid_l) %>%
  mutate(head_shape = head_l/head_w)

# calculate averages across the 10 observations / individual
sperm_means <- sperm_df %>%
  group_by(ID) %>%
  reframe(species = species, avg_head_l = mean(head_l), avg_head_w = mean(head_w), avg_mid_l = mean(mid_l), avg_head_shape = mean(head_shape)) %>%
  distinct()

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
sperm_motility_df <- read_csv("/Users/aguillon/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Xcortezi_Xbirchmanni_reproductive_barriers/Data/Sperm_imaging/Data/mean_motilities.csv", 
                              col_names=TRUE)

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

