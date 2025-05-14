# This script contains the R code to assess ancestry distributions in CHPL
#
# Authors: Aguillon SM, et al.
# Year: 2025
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: Nature Ecology and Evolution
# DOI: https://doi.org/10.1038/s41559-025-02669-9
#
# Edited date: 5 Oct 2023
#
# Please cite the paper if you use these scripts
#


# load packages
library(tidyverse)
library(cowplot)
library(diptest)

# load data
contemporary_HI <- read_tsv("./data/CHPL-contemporary-ancestry.txt", col_names=TRUE)
historic_HI <- read_tsv("./data/CHPL-historic-ancestry.txt", col_names=TRUE)

# plot colors
corcol <- "#359D6F"
bircol <- "#3A6BC6"



##### CONTEMPORARY DATA #####

## Hartigan's dip test for unimodality
dip.test(contemporary_HI$hybrid_index, simulate.p.value=TRUE)
#Hartigans' dip test for unimodality / multimodality with simulated p-value (based
#	on 2000 replicates)
#
#data:  contemporary_HI$hybrid_index
#D = 0.16559, p-value < 2.2e-16
#alternative hypothesis: non-unimodal, i.e., at least bimodal



## summarize ancestry within each "cluster" of individuals
# classify three clusters of individuals: birchmanni, cortezi-like, hybrids
contemporary_HI <- contemporary_HI %>% 
  mutate(cluster=ifelse(hybrid_index>0.6,"cortezi-like",ifelse(hybrid_index<0.1,"birchmanni","hybrid"))) %>%
  group_by(cluster)
# average ancestry in cortezi-like hybrids for each year
summarize(contemporary_HI,
          n=n(),
          prop_total=n/nrow(contemporary_HI),
          avg_ancestry=mean(hybrid_index),
          sd_ancestry=sd(hybrid_index),
          mtDNA_xcor=sum(mtDNA==2),
          mtDNA_xbir=sum(mtDNA==0))
## A tibble: 3 × 7
#cluster          n prop_total avg_ancestry sd_ancestry mtDNA_xcor mtDNA_xbir
#<chr>        <int>      <dbl>        <dbl>       <dbl>      <int>      <int>
#1 birchmanni     189    0.618         0.0189     0.00571          0        189
#2 cortezi-like   115    0.376         0.757      0.0174         115          0
#3 hybrid           2    0.00654       0.244      0.0820           1          1


# histogram of contemporary ancestries (Figure 1B)
contemporary_panel <- ggplot(data=contemporary_HI) + 
  geom_histogram(aes(x=hybrid_index, fill=cut(hybrid_index,5)), bins=30) + 
  scale_fill_manual(values=c(bircol, "black", corcol)) + 
  #xlab("Proportion genome cortezi-derived") + 
  ylab("Count") +
  coord_cartesian(xlim = c(0, 1.0)) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black"))



##### HISTORIC DATA #####

# classify three clusters of individuals: birchmanni, cortezi-like, hybrids
historic_HI <- historic_HI %>% 
  mutate(cluster=ifelse(hybrid_index>0.6,"cortezi-like",ifelse(hybrid_index<0.1,"birchmanni","hybrid"))) %>%
  group_by(year,cluster)
# average ancestry in cortezi-like hybrids for each year
summarize(historic_HI,
          n=n(),
          avg_ancestry=mean(hybrid_index),
          sd_ancestry=sd(hybrid_index),
          mtDNA_xcor=sum(mtDNA==2),
          mtDNA_xbir=sum(mtDNA==0))
## A tibble: 7 × 7
## Groups:   year [3]
#year cluster          n avg_ancestry sd_ancestry mtDNA_xcor mtDNA_xbir
#<dbl> <chr>        <int>        <dbl>       <dbl>      <int>      <int>
#1  2003 birchmanni       4       0.0145     0.00411          0          4
#2  2003 cortezi-like     7       0.782      0.0122           7          0
#3  2006 birchmanni       1       0.0120    NA                0          1
#4  2006 cortezi-like    19       0.796      0.0351          19          0
#5  2006 hybrid           1       0.523     NA                1          0
#6  2017 birchmanni      22       0.0130     0.00513          0         22
#7  2017 cortezi-like    19       0.765      0.0145          19          0

## Hartigan's dip test for unimodality overall
dip.test(historic_HI$hybrid_index, simulate.p.value=TRUE)
#Hartigans' dip test for unimodality / multimodality with simulated p-value (based on
#	2000 replicates)
#
#data:  historic_HI$hybrid_index
#D = 0.17945, p-value < 2.2e-16
#alternative hypothesis: non-unimodal, i.e., at least bimodal



### separate analyses for 2003
dip.test(filter(historic_HI,year==2003)$hybrid_index, simulate.p.value=TRUE)
#Hartigans' dip test for unimodality / multimodality with simulated p-value (based on 2000
#	replicates)
#
#data:  filter(historic_HI, year == 2003)$hybrid_index
#D = 0.17961, p-value = 0.001
#alternative hypothesis: non-unimodal, i.e., at least bimodal

### separate analyses for 2006
dip.test(filter(historic_HI,year==2006)$hybrid_index, simulate.p.value=TRUE)
#Hartigans' dip test for unimodality / multimodality with simulated p-value (based on 2000
#	replicates)
#
#data:  filter(historic_HI, year == 2006)$hybrid_index
#D = 0.042989, p-value = 0.996
#alternative hypothesis: non-unimodal, i.e., at least bimodal

### separate analyses for 2017
dip.test(filter(historic_HI,year==2017)$hybrid_index, simulate.p.value=TRUE)
#Hartigans' dip test for unimodality / multimodality with simulated p-value (based on 2000
#	replicates)
#
#data:  filter(historic_HI, year == 2017)$hybrid_index
#D = 0.21626, p-value < 2.2e-16
#alternative hypothesis: non-unimodal, i.e., at least bimodal





# histogram of historic ancestries (Figure 1C)
historic_panel <- ggplot() + 
  geom_density(data=filter(historic_HI, year==2003), aes(x=hybrid_index), adjust=1/5, color="gray75", linewidth=1) +
  geom_density(data=filter(historic_HI, year==2006), aes(x=hybrid_index), adjust=5, color="gray50", linewidth=1) +
  geom_density(data=filter(historic_HI, year==2017), aes(x=hybrid_index), adjust=1/5, color="black", linewidth=1) + 
  xlab("Proportion genome X. cortezi-derived") + 
  ylab("Density") +
  xlim(0.0,1.0) +
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black"), strip.text.x=element_text(size=10,face="bold"))



#### combine plots for in-text figure ####
plot_grid(contemporary_panel, historic_panel, align="v", nrow=2, rel_heights=c(1,0.75))
# export as 4" x 3.5" PDF
# in Illustrator: combined into larger multi-panel figure with distribution map, added legend and details to each panel, italicized "X. cortezi"
  

#### separate historic plots for supplement #####   
ggplot(data=historic_HI) + 
  geom_histogram(aes(x=hybrid_index, fill=cut(hybrid_index,5)), binwidth=0.05) + 
  facet_grid(cols=vars(year)) + 
  scale_fill_manual(values=c(bircol, "black", corcol)) + 
  xlab("Proportion genome X. cortezi-derived") + 
  ylab("Count") + 
  coord_cartesian(xlim = c(0, 1.0)) +
  theme_bw() + 
  theme(legend.position="none", axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"), strip.text.x=element_text(size=12,face="bold"))
# export as 6.5" x 4" PDF
# in Illustrator: removed excess 0s in x-axis tick labels, italicized "X. cortezi"


