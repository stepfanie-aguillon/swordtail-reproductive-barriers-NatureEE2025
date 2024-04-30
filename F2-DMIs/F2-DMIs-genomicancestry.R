# This script contains the R code to work with DMIs and average ancestry in F2s
#
# Authors: Aguillon SM, et al.
# Year: 2024
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: TBD
# bioRxiv DOI: https://doi.org/10.1101/2024.04.16.589374
#
# Edited date: 8 March 2024
#
# Please cite the paper if you use these scripts
#


# load packages
library(tidyverse)
library(cowplot)
library(ggrepel)

# load data
## selection coefficients from ABC 
ndufs5_sel <- read_tsv("./data/ABC_simulations_accepted_F2s_Feb2024_ch13.txt", col_names=TRUE)
ndufa13_sel <- read_tsv("./data/ABC_simulations_accepted_F2s_Feb2024_ch6.txt" ,col_names=TRUE)
chr7_sel <- read_tsv("./data/ABC_simulations_accepted_F2s_Feb2024_ch14.txt", col_names=TRUE)
chr14_sel <- read_tsv("./data/ABC_simulations_accepted_F2s_Feb2024_ch7.txt" ,col_names=TRUE)

# load data
## Assessment of known hybrid incompatibilities in F2 embryos
actual_F2s <- read_tsv("./F2-sample-IDs.txt", col_names=TRUE)

F2_DMI_df1 <- read_tsv("./data/selected_DMI_XcorXbir-PF2-genotypes-110623.txt_transposed", col_names=TRUE) %>%
  rename(ID = chr, chr13_2112372 = "chr-13", chr6_12520420 = "chr-06", mt_5159 = mitochondria) %>%
  filter(!row_number() %in% c(1)) %>%
  filter(ID %in% actual_F2s$ID)
F2_DMI_df2 <- read_tsv("./data/selected_DMI_XcorXbir_genotypes-S255-021924.txt_transposed", col_names=TRUE) %>%
  rename(ID = chr, chr13_2112372 = "chr-13", chr6_12520420 = "chr-06", mt_5159 = mitochondria) %>%
  filter(!row_number() %in% c(1)) %>%
  filter(ID %in% actual_F2s$ID)

F2_DMI_df <- rbind(F2_DMI_df1, F2_DMI_df2)


# data structure
head(F2_DMI_df)
## A tibble: 6 × 4
#ID                                             chr13_2112372 chr6_12520420 mt_5159
#<chr>                                                  <dbl>         <dbl>   <dbl>
#1 XcorXbir-F1-5-IX-23-S170-L-GGY-E01-PI.R1.fastq             0             1       2
#2 XcorXbir-F1-5-IX-23-S170-L-GGY-E02-PI.R1.fastq             0             1       2
#3 XcorXbir-F1-5-IX-23-S170-L-GGY-E03-PI.R1.fastq             0             1       2

# about the 0, 1, 2 values
#0 = Xbir homozygote (0 copies Xcor)
#1 = heterozygote (1 copy Xcor)
#2 = Xcor homozygote (2 copies Xcor)


# filter just to the embryos, indicate ones that are potentially incompatible based on phenotype
embryos <- F2_DMI_df %>%
  filter(str_detect(ID, "-E")) %>%
  mutate(DMI = ifelse(str_detect(ID, "-PI"), "Y", "N"))



# chr 13, gene = ndufs5
DMI_summary1 <- embryos %>%
  group_by(DMI, chr13_2112372) %>%
  summarize(count = n()) %>% 
  mutate(prop = ifelse(DMI=="N", count/(64+29), count/33))
DMI_summary1
## A tibble: 3 × 4
## Groups:   DMI [2]
#DMI   chr13_2112372 count  prop
#<chr>         <dbl> <int> <dbl>
#1 N                 1    64 0.688
#2 N                 2    29 0.312
#3 Y                 0    33 1   


# chi-squared test for DMI = N
chisq.test(c(0,64,29), p=c(0.25,0.5,0.25))
#Chi-squared test for given probabilities
#
#data:  c(0, 64, 29)
#X-squared = 31.258, df = 2, p-value = 1.631e-07

# expected proportions:
(64+29)*0.25 #23.25
(64+29)*0.5 #46.5

# chi-squared test for DMI = Y
chisq.test(c(33,0,0), p=c(0.25,0.5,0.25))
#Chi-squared test for given probabilities
#
#data:  c(33, 0, 0)
#X-squared = 99, df = 2, p-value < 2.2e-16
chisq.test(c(33,0,0), p=c(0.25,0.5,0.25))$p.value
#3.179971e-22

# expected proportions:
33*0.25 #8.25
33*0.5 #16.5


# chr 6, gene = ndufa13
DMI_summary2 <- embryos %>%
  group_by(DMI, chr6_12520420) %>%
  summarize(count = n()) %>% 
  mutate(prop = ifelse(DMI=="N", count/(26+51+16), count/(9+17+7)))
DMI_summary2
## A tibble: 6 × 4
## Groups:   DMI [2]
#DMI   chr6_12520420 count  prop
#<chr>         <dbl> <int> <dbl>
#1 N                 0    26 0.280
#2 N                 1    51 0.548
#3 N                 2    16 0.172
#4 Y                 0     9 0.273
#5 Y                 1    17 0.515
#6 Y                 2     7 0.212

# chi-squared test for DMI = Y
chisq.test(c(26,51,16), p=c(0.25,0.5,0.25))
#Chi-squared test for given probabilities
#
#data:  c(26, 51, 16)
#X-squared = 3.0215, df = 2, p-value = 0.2207

# chi-squared test for DMI = N
chisq.test(c(9,17,7), p=c(0.25,0.5,0.25))
#Chi-squared test for given probabilities
#
#data:  c(9, 17, 7)
#X-squared = 0.27273, df = 2, p-value = 0.8725





##### analysis in adult samples
adults <- F2_DMI_df %>%
  filter(!str_detect(ID, "-E"))


# chr 13, gene = ndufs5
DMI_summary1_adults <- adults %>%
  group_by(chr13_2112372) %>%
  summarize(count_chr13 = n(), prop_chr13=count_chr13/(113+50)) 
DMI_summary1_adults
## A tibble: 2 × 3
#chr13_2112372 count_chr13 prop_chr13
#<dbl>       <int>      <dbl>
#1             1         113      0.693
#2             2          50      0.307

# chi-squared test
chisq.test(c(0,113,50), p=c(0.25,0.5,0.25))
#Chi-squared test for given probabilities
#
#data:  c(0, 113, 50)
#X-squared = 55.025, df = 2, p-value = 1.126e-12


# chr 6, gene = ndufa13
DMI_summary2_adults <- adults %>%
  group_by(chr6_12520420) %>%
  summarize(count_chr6 = n(), prop_chr6=count_chr6/(23+95+45)) 
DMI_summary2_adults
## A tibble: 3 × 3
#chr6_12520420 count_chr6 prop_chr6
#<dbl>      <int>     <dbl>
#1             0         23     0.141
#2             1         95     0.583
#3             2         45     0.276

# chi-squared test
chisq.test(c(23,95,45), p=c(0.25,0.5,0.25))
#Chi-squared test for given probabilities
#
#data:  c(23, 95, 45)
#X-squared = 10.411, df = 2, p-value = 0.005486


# combine summary dfs for two genes
# need to rename to make it work
adults_1 <- DMI_summary1_adults %>%
  rename(Xcor = "chr13_2112372") %>%
  select(Xcor, count_chr13)
adults_2 <- DMI_summary2_adults %>%
  rename(Xcor = "chr6_12520420") %>%
  select(Xcor, count_chr6)
combined_df <- full_join(adults_1, adults_2,by="Xcor") %>%
  pivot_longer(!Xcor, names_to = "gene_ID", values_to = "count") %>% 
  mutate(prop = ifelse(gene_ID == "count_chr13", count/(113+50), count/(23+95+45)))


# assess individuals that are homozygous Xbir for ndufs5 across all F2 samples
ndufs5 <- F2_DMI_df %>%
  filter(chr13_2112372 == 0) %>%
  mutate(age = ifelse(str_detect(ID, "-E"), "embryo", "adult"))
nrow(ndufs5)
# 33 samples, ALL are the identified embryos that were IDed as potentially incompatible based on the phenotype

# assess individuals that are homozygous Xbir for ndufa13 across all F2 samples
ndufa13 <- F2_DMI_df %>%
  filter(chr6_12520420 == 0) %>%
  mutate(age = ifelse(str_detect(ID, "-E"), "embryo", "adult"))
nrow(ndufa13)
# 58 samples total, only 23 are from sequenced adults, remaining 35 are from embryos dissected/sequenced before birth
35/58 #0.6034483
23/58 #0.3965517




##### figures for ndufs5 and ndufa13 for main text and supplement #####

# function to determine maximum a posteriori estimate for ABC simulations
posterior.mode = function(x){
  if(!is.na(x[1])){ 
    x.max = max(x)
    x.min = min(x)
    dres <- density( x, from = x.min, to = x.max)
    modeParam <- dres$x[which.max(dres$y)]
  } else if(is.na(x[1])){
    modeParam <- NA
  }
  return(modeParam)
}
#to run: posterior.mode(df$variable)


# ndufs5 color = #E6AB02
# ndufa13 color = #7570B3


### main text figure combining ndufs5 in embryos and ndufa13 in adults
ndufs5_ndufa13_plot <- ggplot() + 
  geom_count(data=DMI_summary1, aes(x=DMI, y=as.character(chr13_2112372), size=count), color="#E6AB02") + 
  geom_label_repel(data=DMI_summary1, aes(x=DMI, y=as.character(chr13_2112372), label=paste0(round(prop*100,1),"%")), nudge_y = -0.3, segment.color="NA", size=2, label.padding=0.1) +
  geom_count(data=filter(combined_df, gene_ID=="count_chr6"), aes(y=as.character(Xcor), x=gene_ID, size=count), color="#7570B3") + 
  geom_label_repel(data=filter(combined_df, gene_ID=="count_chr6"), aes(x=gene_ID, y=as.character(Xcor), label=paste0(round(prop*100,1),"%")), nudge_y = -0.3, segment.color="NA", size=2, label.padding=0.1) +
  scale_size_area(limits=c(0,120)) +
  ylab("Genotype") +
  xlab("update axis!") + 
  scale_y_discrete(labels=c("BB", "BC", "CC")) +
  scale_x_discrete(limits=c("Y", "N", "count_chr6"), labels=c("stalled embryo", "normal embryo", "adult")) +
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black"))
# to be combined with ancestry plots

### main text figure combining ABC results for ndufs5 and ndufa13
ndufs5_ndufa13_sel_plot <- ggplot() +
  geom_density(data=ndufs5_sel, aes(x=s, y=after_stat(density*nrow(ndufs5_sel))/7500), fill="#E6AB02", color=NA) + 
  geom_vline(xintercept=posterior.mode(ndufs5_sel$s), color="black", linetype="dashed") +
  geom_density(data=ndufa13_sel, aes(x=s, y=after_stat(density*nrow(ndufa13_sel))/7500), fill="#7570B3", color=NA) +
  geom_vline(xintercept=posterior.mode(ndufa13_sel$s), color="black", linetype="dashed") +
  xlim(0,1.0) + 
  ylab("Density") +
  xlab("Selection coefficient") +
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black"))
# to be combined with ancestry plots

posterior.mode(ndufs5_sel$s) #0.9953891
round(quantile(ndufs5_sel$s,c(0.025,0.975)),3)
#2.5%     97.5% 
#0.933 1.000 

posterior.mode(ndufa13_sel$s) #0.5314902
round(quantile(ndufa13_sel$s,c(0.025,0.975)),2)
#2.5%     97.5% 
#0.201 0.694 



### supplemental figure for ndufa13 in embryos
ndufa13_plot <- ggplot(data=DMI_summary2) + 
  geom_count(aes(x=DMI, y=as.character(chr6_12520420), size=count), color="#7570B3") + 
  geom_label_repel(aes(x=DMI, y=as.character(chr6_12520420), label=paste0(round(prop*100,1),"%")), nudge_y = -0.3, segment.color="NA", size=3) +
  scale_size_area(limits=c(0,51)) +
  ylab("Genotype") +
  xlab("Embryonic development") + 
  scale_x_discrete(limits=c("Y", "N"), labels=c("stalled", "normal")) +
  scale_y_discrete(labels=c("BB", "BC", "CC")) +
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"))
# export as 3" x 3" PDF






##### average ancestry across chromosomes

# load data
avg_ancestry <- read_tsv("./data/XcorXbir-F2-fullgenotypes.txt_transposed_ancestry_by_site", col_names=TRUE) %>%
  mutate(sig = hybrid_index > 0.58)

# average genome-wide
mean(avg_ancestry$hybrid_index) #0.5063693


# Xcor ancestry significance threshold from simulations 
HI_min = 0.42 
HI_max = 0.58
#more conservative 1%
#HI_min = 0.40
#HI_max = 0.60

# regions that have more or less ancestry than the simulated significance threshold
more_xcor <- avg_ancestry %>%
  filter(hybrid_index >= HI_max)

less_xcor <- avg_ancestry %>%
  filter(hybrid_index <= HI_min)



# ndufs5 color = #E6AB02
# ndufa13 color = #7570B3

# chr 13, ndufs5
#### plot the results for multipanel figure ####
chr13_ancestry <- avg_ancestry %>%
  filter(group=="chr-13") %>%
  filter(!between(position, 20538312, 20875924))

chr13_ancestry_plot <- ggplot(data=chr13_ancestry) + 
  geom_ribbon(aes(x=position/1e6, ymin=HI_min, ymax=HI_max), fill="light gray") +
  geom_point(aes(x=position/1e6, y=hybrid_index), color="black", size=0.5) + 
  geom_point(aes(x=2112372/1e6, y=0.3), color="#E6AB02", shape=17, size=5) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  xlab("Chromosome 13 position (Mb)") +
  ylab("Proportion genome X. cortezi-derived") + 
  ylim(0.3,0.7) + 
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black"))


# chr 6, ndufa13
#### plot the results for multipanel figure ####
chr6_ancestry <- avg_ancestry %>%
  filter(group=="chr-06") %>%
  filter(!between(position, 7754852, 8085477)) %>%
  filter(!between(position, 23899215, 24123182))

chr6_ancestry_plot <- ggplot(data=chr6_ancestry) + 
  geom_ribbon(aes(x=position/1e6, ymin=HI_min, ymax=HI_max), fill="light gray") +
  geom_point(aes(x=position/1e6, y=hybrid_index), size=0.5, color="black") + 
  geom_point(aes(x=12520420/1e6, y=0.3), color="#7570B3", shape=17, size=5) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  xlab("Chromosome 6 position (Mb)") +
  ylab("Proportion genome X. cortezi-derived") + 
  ylim(0.3,0.7) + 
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black"))



#### PLOT MULTIPANEL FOR FIGURE 4 ####
plot_grid(ndufs5_ndufa13_plot, chr13_ancestry_plot, ndufs5_ndufa13_sel_plot, chr6_ancestry_plot, nrow=2, 
          rel_widths=c(1, 1.75, 1, 1.75))
# export as 4.5" tall x 7" wide
# in Illustrator: fix x-axis for first panel, add selection coefficients to density plots, add gene labels to ancestry plots
# italicize species names, fix y-axis for ancestry plots, combine with other figures, add panel labels






#### chr 1 for comparison figure ####
ggplot(data=filter(avg_ancestry, group=="chr-01")) + 
  geom_ribbon(aes(x=position/1e6, ymin=HI_min, ymax=HI_max), fill="light gray") +
  geom_point(aes(x=position/1e6, y=hybrid_index), color="black", size=0.5) + 
  geom_hline(yintercept=0.5, linetype="dashed") +
  xlab("Chromosome 1 position (Mb)") +
  ylab("Proportion genome X. cortezi-derived") + 
  ylim(0.3,0.70) + 
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black"))
# export as 4" high x 5" wide
# in Illustrator: italicize species name




#### PLOT MULTIPANEL FOR FIGURE 5 - additional segregation distorters ####


#### chr 7 region
chr7_ancestry <- avg_ancestry %>%
  filter(group=="chr-07") %>%
  filter(!between(position, 24792083, 25737000)) 

chr7_ancestry_plot <- ggplot() + 
  geom_ribbon(data=chr7_ancestry, aes(x=position/1e6, ymin=HI_min, ymax=HI_max), fill="light gray") +
  geom_point(data=chr7_ancestry, aes(x=position/1e6, y=hybrid_index, color=sig), size=0.5) + 
  scale_color_manual(values=c("black", "red")) +
  geom_point(data=filter(chr7_ancestry, sig=="TRUE"), aes(x=position/1e6, y=hybrid_index), size=0.75, color="red") +
  geom_point(aes(x=22576066/1e6, y=0.385), color="red", shape=17, size=3) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  geom_segment(aes(x=16781454/1e6, xend=28711584/1e6, y=0.385, yend=0.385), color="red", linewidth=1) +
  xlab("Chromosome 7 position (Mb)") +
  ylab("Proportion genome X. cortezi-derived") + 
  scale_y_continuous(limits=c(0.375,0.625), breaks=c(0.4,0.5,0.6)) +
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black"))

chr7_sel_plot <- ggplot(data=chr7_sel) + 
  geom_density(aes(x=s), fill="dark gray", color=NA) + 
  geom_vline(xintercept=posterior.mode(chr7_sel$s), color="black", linetype="dashed") +
  xlim(0,1.0) + 
  ylim(0,4) +
  ylab("Density") +
  xlab("Selection coefficient") +
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black"))

posterior.mode(chr7_sel$s) #0.4689934
round(quantile(chr7_sel$s,c(0.025,0.975)),3)
#2.5%     97.5% 
#0.177 0.626 


#### chr 14 region
chr14_ancestry <- avg_ancestry %>%
  filter(group=="chr-14") %>%
  filter(!between(position, 18269172, 18539743))

chr14_ancestry_plot <- ggplot() + 
  geom_ribbon(data=chr14_ancestry, aes(x=position/1e6, ymin=HI_min, ymax=HI_max), fill="light gray") +
  geom_point(data=chr14_ancestry, aes(x=position/1e6, y=hybrid_index, color=sig), size=0.5) + 
  scale_color_manual(values=c("black", "red")) +
  geom_point(data=filter(chr14_ancestry, sig=="TRUE"), aes(x=position/1e6, y=hybrid_index), size=0.75, color="red") +
  geom_point(aes(x=9410308/1e6, y=0.385), color="red", shape=17, size=3) +
  geom_segment(aes(x=6569477/1e6, xend=11511130/1e6, y=0.385, yend=0.385), color="red", linewidth=1) + 
  geom_hline(yintercept=0.5, linetype="dashed") +
  xlab("Chromosome 14 position (Mb)") +
  ylab("Proportion genome X. cortezi-derived") + 
  scale_y_continuous(limits=c(0.375,0.625), breaks=c(0.4,0.5,0.6)) +
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black"))

chr14_sel_plot <- ggplot(data=chr14_sel) + 
  geom_density(aes(x=s), fill="dark gray", color=NA) + 
  geom_vline(xintercept=posterior.mode(chr14_sel$s), color="black", linetype="dashed") +
  xlim(0,1.0) + 
  ylim(0,4) +
  ylab("Density") +
  xlab("Selection coefficient") +
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black"))

posterior.mode(chr14_sel$s) #0.5269006
round(quantile(chr14_sel$s,c(0.025,0.975)),3)
#2.5%     97.5% 
#0.229 0.703 


#### PLOT MULTIPANEL FOR FIGURE 5 ####
plot_grid(chr7_ancestry_plot, chr7_sel_plot, chr14_ancestry_plot, chr14_sel_plot, nrow=2, 
          rel_widths=c(1, 1, 1, 1))
# export as 3.5" tall x 3.5" wide
# in Illustrator: fix y-axis for ancestry panels, add panel labels, italicize species names, add selection coefficients








### plotting deserts and islands from CHPL and STAC

# load data
minor_parent_data <- read_tsv("./data/xbir_pacbio2023_allChrs_10kb_windows_recRate_codingBP_conservedBP_everything_minPar.txt", col_names=TRUE)

# load average ancestry to plot significant regions on top of plot
avg_ancestry <- read_tsv("./data/XcorXbir-F2-fullgenotypes.txt_transposed_ancestry_by_site", col_names=TRUE) %>%
  mutate(sig = hybrid_index > 0.58)
chr7_sig <- avg_ancestry %>%
  filter(group=="chr-07") %>%
  filter(!between(position, 24792083, 25737000)) %>%
  filter(sig == "TRUE")
chr14_sig <- avg_ancestry %>%
  filter(group=="chr-14") %>%
  filter(!between(position, 18269172, 18539743)) %>%
  filter(sig == "TRUE")

## chr 7 plot
chr7_plot <- ggplot() + 
  geom_point(data=filter(minor_parent_data, chrID=="chr-07"), aes(x=start/1e6, y=1-CHPL2021), color="black", size=0.5, alpha=0.5) + 
  geom_point(data=filter(minor_parent_data, chrID=="chr-07"), aes(x=start/1e6, y=1-STAC2020), color="gray", size=0.5, alpha=0.5) + 
  geom_point(data=chr7_sig, aes(x=position/1e6, y=1.05), color="red", shape=15, size=0.5) +
  geom_segment(aes(x=16781454/1e6, xend=28711584/1e6, y=1.05, yend=1.05), color="red", linewidth=0.25) +
  xlab("Chromosome 7 position (Mb)") +
  ylab("Proportion genome X. cortezi-derived") + 
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"))

## chr 14 plot
chr14_plot <- ggplot() + 
  geom_point(data=filter(minor_parent_data, chrID=="chr-14"), aes(x=start/1e6, y=1-CHPL2021), color="black", size=0.5, alpha=0.5) + 
  geom_point(data=filter(minor_parent_data, chrID=="chr-14"), aes(x=start/1e6, y=1-STAC2020), color="gray", size=0.5, alpha=0.5) + 
  geom_point(data=chr14_sig, aes(x=position/1e6, y=1.05), color="red", shape=15, size=0.5) +
  geom_segment(aes(x=6569477/1e6, xend=11511130/1e6, y=1.05, yend=1.05), color="red", linewidth=0.25) + 
  xlab("Chromosome 14 position (Mb)") +
  ylab("Proportion genome X. cortezi-derived") + 
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"))

#### multipanel for supplement ####
plot_grid(chr7_plot, chr14_plot, nrow=2)
# export as 6" tall x 6" wide
# in Illustrator: italicize species names, fix y-axis, add panel labels









#### plotting h for the ABC simulations ####

panel1 <- ggplot(data=ndufs5_sel) + 
  geom_density(aes(x=h, y=after_stat(density*nrow(ndufs5_sel))/1200), fill="#E6AB02", color=NA) + 
  geom_vline(xintercept=posterior.mode(ndufs5_sel$h), color="black", linetype="dashed") +
  xlim(0,1.0) + 
  ylab("Density") +
  xlab("Dominance coefficient") +
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black")) +
  theme(axis.title.x=element_blank())

posterior.mode(ndufs5_sel$h) #0.02748369
round(quantile(ndufs5_sel$h,c(0.025,0.975)),3)
#2.5%       97.5% 
#0.004 0.267


panel2 <- ggplot(data=ndufa13_sel) + 
  geom_density(aes(x=h, y=after_stat(density*nrow(ndufa13_sel))/8000), fill="#7570B3", color=NA) + 
  geom_vline(xintercept=posterior.mode(ndufa13_sel$h), color="black", linetype="dashed") +
  xlim(0,1.0) + 
  ylab("Density") +
  xlab("Dominance coefficient") +
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black")) +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank())

posterior.mode(ndufa13_sel$h) #0.04892079
round(quantile(ndufa13_sel$h,c(0.025,0.975)),2)
#2.5%       97.5% 
#0.008 0.606  

panel3 <- ggplot(data=chr7_sel) + 
  geom_density(aes(x=h, y=after_stat(density*nrow(chr7_sel))/2200), fill="dark gray", color=NA) + 
  geom_vline(xintercept=posterior.mode(chr7_sel$h), color="black", linetype="dashed") +
  xlim(0,1.0) + 
  ylab("Density") +
  xlab("Dominance coefficient") +
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black"))

posterior.mode(chr7_sel$h) #0.9410835
round(quantile(chr7_sel$h,c(0.025,0.975)),3)
#0.398 0.994 

panel4 <- ggplot(data=chr14_sel) + 
  geom_density(aes(x=h, y=after_stat(density*nrow(chr14_sel))/1700), fill="dark gray", color=NA) + 
  geom_vline(xintercept=posterior.mode(chr14_sel$h), color="black", linetype="dashed") +
  xlim(0,1.0) + 
  ylab("Density") +
  xlab("Dominance coefficient") +
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black")) +
  theme(axis.title.y=element_blank())

posterior.mode(chr14_sel$h) #0.5182547
round(quantile(chr14_sel$h,c(0.025,0.975)),3)
#2.5%      97.5% 
#0.055 0.892 

plot_grid(panel1, panel2, panel3, panel4, nrow=2)
# export as 4" tall x 6" wide
# in Illustrator: add panel labels, add h values (similar to s values)








##### Genes in regions on chromosome 7 and chromosome 14

# .gff file header info
gff_header <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

#load gff to look at genes in segregation distorted regions
xbir_gff <- read_tsv("./data/xbir-COAC-16-VIII-22-M_v2023.1.gff", col_names=FALSE) %>%
  rename_with(~gff_header)

# region on chromosome 7
# low = 16781454, high = 30564405
chr7_genes <- xbir_gff %>%
  filter(seqname == "chr-07" & feature =="gene") %>%
  filter(end <= 30564405) %>%
  filter(start >= 16781454)
#write.table(chr7_genes, "./chr7-segdis-genelist.txt", quote=FALSE, sep="\t", row.names=FALSE)

# region on chromosome 14
# low = 6569477, high = 11511130
chr14_genes <- xbir_gff %>%
  filter(seqname == "chr-14" & feature =="gene") %>%
  filter(end <= 11511130) %>%
  filter(start >= 6569477)
#write.table(chr14_genes, "./chr14-segdis-genelist.txt", quote=FALSE, sep="\t", row.names=FALSE)
