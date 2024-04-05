# This script contains the R code to assess ancestry of cross-cluster hybrids
#
# Authors: Aguillon SM, et al.
# Year: 2024
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: TBD
# bioRxiv DOI: TBD
#
# Edited date: 5 April 2024
#
# Please cite the paper if you use these scripts
#


# load packages
library(tidyverse)
library(cowplot)


##### ANCESTRY TRACTS FOR CROSS-CLUSTER HYBRIDS #####

# load data
df1_contemporary <- read_tsv("./data/for-ancestry-tracts/CHPL2021_2023v1pacbio_xcorxbir_genotypes.txt_transposed", col_names=TRUE)
df2_contemporary <- read_tsv("./data/for-ancestry-tracts/CHPL-III-22_CHPL-VI-22_CHPL-V-22_XcorXbir-PF1_genotypes.txt_transposed", col_names=TRUE)
df_2006 <- read_tsv("./data/for-ancestry-tracts/CHPL2006_2023v1pacbio_xcorxbir_genotypes.txt_transposed", col_names=TRUE)


### CHROMOSOME 1 PLOT

#sample: swt-CHPL-VI-2006-M04.R1.fastq
df_2006_filter <- df_2006 %>%
  filter(chr=="chr-01") %>%
  select(chr,pos,"swt-CHPL-VI-2006-M04.R1.fastq") %>%
  drop_na()
plot(main="CHPL-VI-2006-M04", x=df_2006_filter$pos/1e6, y=df_2006_filter$`swt-CHPL-VI-2006-M04.R1.fastq`, type="l", lwd=2, col="black", ylim=c(0,2), yaxt="n", xlab="", ylab="# cortezi alleles")
axis(2, at=c(0,1,2))

#sample: CHPL-3-VI-21-68.R1.fastq
df1_contemporary_filter2 <- df1_contemporary %>%
  filter(chr=="chr-01") %>%
  select(chr,pos,"CHPL-3-VI-21-68.R1.fastq") %>%
  drop_na()
plot(main="CHPL-3-VI-21-68", x=df1_contemporary_filter2$pos/1e6, y=df1_contemporary_filter2$`CHPL-3-VI-21-68.R1.fastq`, type="l", lwd=2, col="black", ylim=c(0,2), yaxt="n", xlab="", ylab="# cortezi alleles")
axis(2, at=c(0,1,2))

#sample: CHPL-VI-22-32-R-VOV-M-REDO.R1.fastq
df2_contemporary_filter <- df2_contemporary %>%
  filter(chr=="chr-01") %>%
  select(chr,pos,"CHPL-VI-22-32-R-VOV-M-REDO.R1.fastq") %>%
  drop_na()
plot(main="CHPL-VI-22-32-R-VOV-M", x=df2_contemporary_filter$pos/1e6, y=df2_contemporary_filter$`CHPL-VI-22-32-R-VOV-M-REDO.R1.fastq`, type="l", lwd=2, col="black", ylim=c(0,2), yaxt="n", xlab="Chromosome 1 Position (Mb)", ylab="# cortezi alleles")
axis(2, at=c(0,1,2))


### CHROMOSOME 2 PLOT

#sample: swt-CHPL-VI-2006-M04.R1.fastq
df_2006_filter <- df_2006 %>%
  filter(chr=="chr-02") %>%
  select(chr,pos,"swt-CHPL-VI-2006-M04.R1.fastq") %>%
  drop_na()
plot(main="CHPL-VI-2006-M04", x=df_2006_filter$pos/1e6, y=df_2006_filter$`swt-CHPL-VI-2006-M04.R1.fastq`, type="l", lwd=2, col="black", ylim=c(0,2), yaxt="n", xlab="", ylab="# cortezi alleles")
axis(2, at=c(0,1,2))

#sample: CHPL-3-VI-21-68.R1.fastq
df1_contemporary_filter2 <- df1_contemporary %>%
  filter(chr=="chr-02") %>%
  select(chr,pos,"CHPL-3-VI-21-68.R1.fastq") %>%
  drop_na()
plot(main="CHPL-3-VI-21-68", x=df1_contemporary_filter2$pos/1e6, y=df1_contemporary_filter2$`CHPL-3-VI-21-68.R1.fastq`, type="l", lwd=2, col="black", ylim=c(0,2), yaxt="n", xlab="", ylab="# cortezi alleles")
axis(2, at=c(0,1,2))

#sample: CHPL-VI-22-32-R-VOV-M-REDO.R1.fastq
df2_contemporary_filter <- df2_contemporary %>%
  filter(chr=="chr-02") %>%
  select(chr,pos,"CHPL-VI-22-32-R-VOV-M-REDO.R1.fastq") %>%
  drop_na()
plot(main="CHPL-VI-22-32-R-VOV-M", x=df2_contemporary_filter$pos/1e6, y=df2_contemporary_filter$`CHPL-VI-22-32-R-VOV-M-REDO.R1.fastq`, type="l", lwd=2, col="black", ylim=c(0,2), yaxt="n", xlab="Chromosome 2 Position (Mb)", ylab="# cortezi alleles")
axis(2, at=c(0,1,2))

# export individual plots as 3.25" x 3.25" PDFs
# in Illustrator: combined into a single multi-panel, italicized species names, fixed plot spacing






###### SIMULATIONS OF CROSS-CLUSTER ANCESTRY ######

# plot colors
corcol <- "#359D6F"
bircol <- "#3A6BC6"

# load data from observed contemporary sampling
contemporary_HI <- read_tsv("./data/CHPL-contemporary-ancestry.txt", col_names=TRUE) %>%
  mutate(cluster=ifelse(hybrid_index>0.6,"cortezi-like",ifelse(hybrid_index<0.1,"birchmanni","hybrid")))


# summary from the contemporary data from `CHPL-ancestry-distributions.R`
#details on the different CHPL clusters
## A tibble: 3 × 7
#cluster          n prop_total avg_ancestry sd_ancestry mtDNA_xcor mtDNA_xbir
#<chr>        <int>      <dbl>        <dbl>       <dbl>      <int>      <int>
#1 birchmanni     189    0.618         0.0189     0.00571          0        189
#2 cortezi-like   115    0.376         0.757      0.0174         115          0
#3 hybrid           2    0.00654       0.244      0.0820           1          1



### CORTEZI-LIKE CLUSTER
# mean and variance parameters
piC <- 0.757
pi_varC <- 0.0174^2

# beta parameters, method of moments
aC = piC*((piC*(1-piC))/pi_varC-1)
bC = (1-piC)*((piC*(1-piC))/pi_varC-1)
cortezi_like <- rbeta(1000, shape1 = aC, shape2 = bC) 



### BIRCHMANNI CLUSTER
# mean and variance parameters
piB <- 0.0189
pi_varB <- 0.00571^2

# beta parameters, method of moments
aB = piB*((piB*(1-piB))/pi_varB-1)
bB = (1-piB)*((piB*(1-piB))/pi_varB-1)
birchmanni <- rbeta(1000, shape1 = aB, shape2 = bB)


###### determining if the beta distribution nicely fits the actual data
cor_plotting <- as.data.frame(cortezi_like)
bir_plotting <- as.data.frame(birchmanni)
ggplot() + 
  #xlim(0,1) +
  geom_histogram(data=contemporary_HI, aes(x=hybrid_index, fill=cut(hybrid_index,5)), binwidth=0.025) + 
  geom_density(data=cor_plotting, aes(x=cortezi_like), fill = "gray", alpha=0.5) + 
  geom_density(data=bir_plotting, aes(x=birchmanni), fill = "gray", alpha=0.5) + 
  scale_fill_manual(values=c(bircol, "black", corcol)) + 
  coord_cartesian(xlim=c(0,1.0)) +
  xlab("Proportion genome cortezi-derived") +
  ylab("Count") +
  theme_bw() +
  theme(legend.position="none", axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"))







### Simulating a distribution of potential inter-cluster hybrids 
# averaging values from the beta distributions of the birchmanni cluster and the cortezi-like cluster
# adding in some noise from the known sd between embryos from a female

ancestry_sd <- 0.0047 #from mother/embryo sequencing
df <- c()

for(i in 1:length(cortezi_like)){
  
  mu <- (cortezi_like[i]+birchmanni[i])/2
  df[i] <- rnorm(1,mean=mu,sd=ancestry_sd)
  
}

mean_df <- mean(df) 
#0.3875852
sd_df <- sd(df)
#0.01064775
#two SD above, two SD below
#0.3662897 to 0.4088807

#### plot the full results 
df_plotting <- as.data.frame(df)

ggplot(data=df_plotting) + 
  xlim(0,1) +
  geom_density(aes(x=df)) +
  xlab("Proportion genome cortezi-derived") +
  ylab("Simulated Count") +
  theme_bw() +
  theme(legend.position="none", axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"))




### Simulating a distribution of potential inter-cluster hybrids backcrossing into birchmanni
# averaging values from the beta distributions of the inter-cluster hybrid distribution and the birchmanni cluster 
df2 <- c()

for(i in 1:length(cortezi_like)){
  
  mu2 <- (birchmanni[i]+df[i])/2
  df2[i] <- rnorm(1,mean=mu2,sd=ancestry_sd)
  
}

mean_df2 <- mean(df2) 
#0.2030271
sd_df2 <- sd(df2)
#0.008361953
#two SD above, two SD below
#0.1863032 to 0.219751

#### plot the results
df2_plotting <- as.data.frame(df2)

ggplot(data=df2_plotting) + 
  xlim(0,1) +
  geom_density(aes(x=df2)) +
  xlab("Proportion genome cortezi-derived") +
  ylab("Simulated Count") +
  theme_bw() +
  theme(legend.position="none", axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"))




### Simulating a distribution of potential inter-cluster hybrids backcrossing into cortezi-like cluster
# averaging values from the beta distributions of the inter-cluster hybrid distribution and the cortezi-like cluster 
df3 <- c()

for(i in 1:length(cortezi_like)){
  
  mu3 <- (cortezi_like[i]+df[i])/2
  df3[i] <- rnorm(1,mean=mu3,sd=ancestry_sd)
  
}

mean_df3 <- mean(df3) 
#0.5723761
sd_df3 <- sd(df3)
#0.01432684
#two SD above, two SD below
#0.5437224 to 0.6010298

#### plot the results
df3_plotting <- as.data.frame(df3)

ggplot(data=df3_plotting) + 
  xlim(0,1) +
  geom_density(aes(x=df3)) +
  xlab("Proportion genome cortezi-derived") +
  ylab("Simulated Count") +
  theme_bw() +
  theme(legend.position="none", axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"))




### Simulating a distribution of potential third generation cross between an inter-cluster hybrid and a backcross into the birchmanni cluster
# averaging values from the beta distributions of the inter-cluster hybrid distribution and the backross with the birchmanni cluster 
df4 <- c()

for(i in 1:length(df2)){
  
  mu4 <- (df2[i]+df[i])/2
  df4[i] <- rnorm(1,mean=mu4,sd=ancestry_sd)
  
}

mean_df4 <- mean(df4) 
#0.2964181
sd_df4 <- sd(df4)
#0.009839832
#two SD above, two SD below
#0.2767384 to 0.3160978

#### plot the results
df4_plotting <- as.data.frame(df4)

ggplot(data=df4_plotting) + 
  xlim(0,1) +
  geom_density(aes(x=df4)) +
  xlab("Proportion genome cortezi-derived") +
  ylab("Simulated Count") +
  theme_bw() +
  theme(legend.position="none", axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"))


# save simulation results for future use
write.table(df, "./inter-cluster-cross-sim.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(df2, "./birchmanni-backcross-sim.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(df3, "./cortezi-like-backcross-sim.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(df4, "./third-gen-inter-cross-sim.txt", quote=FALSE, sep="\t", row.names=FALSE)




#### plot the combined results for the supplement ####

# load simulation results for plotting
df_plotting <- as.data.frame(read_tsv("inter-cluster-cross-sim.txt", col_names=TRUE))
df2_plotting <- as.data.frame(read_tsv("birchmanni-backcross-sim.txt", col_names=TRUE))
df3_plotting <- as.data.frame(read_tsv("cortezi-like-backcross-sim.txt", col_names=TRUE))
df4_plotting <- as.data.frame(read_tsv("third-gen-inter-cross-sim.txt", col_names=TRUE))


### supplemental figure showing observed data (as histograms) and simulation results (as density plots)
ggplot() + 
  #simulations plotted with geom_density are scaled to fit more in view with the observed data
  geom_density(data=df_plotting, aes(x=df), fill="black", color="black") +
  geom_density(data=df2_plotting, aes(x=df2), fill="#c3d2ee", color="#c3d2ee") +
  geom_density(data=df3_plotting, aes(x=df3), fill="#c0ead7", color="#c0ead7") +
  geom_density(data=df4_plotting, aes(x=df4), fill="light gray", color="light gray") +
  #observed data plotted with geom_histogram (as in main text)
  geom_histogram(data=contemporary_HI, aes(x=hybrid_index, fill=cut(hybrid_index,5)), binwidth=0.025) + 
  scale_fill_manual(values=c(bircol, "red", corcol)) + 
  xlab("Proportion genome cortezi-derived") +
  ylab("Count") +
  coord_cartesian(xlim = c(0,1.0)) +
  theme_bw() +
  theme(legend.position="none", axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"))
# export as 5" x 4" PDF
# in Illustrator: italicized species names, added legend and details to each panel 







##################

# just a short analysis to see how many individuals of cross-cluster ancestry
# we would expect to find in a sampling of the size of our mother-embryo sequencing

# load data and summarize
contemporary_HI <- read_tsv("./data/CHPL-contemporary-ancestry.txt", col_names=TRUE) %>% 
  mutate(cluster=ifelse(hybrid_index>0.6,"cortezi-like",ifelse(hybrid_index<0.1,"birchmanni","hybrid"))) %>%
  group_by(cluster)
# average ancestry in cortezi-like hybrids for each year
summarize(contemporary_HI,
          n=n(),
          avg_ancestry=mean(hybrid_index))
## A tibble: 3 × 3
#cluster          n avg_ancestry
#<chr>        <int>        <dbl>
#1 birchmanni     189       0.0189
#2 cortezi-like   115       0.757 
#3 hybrid           2       0.244 


# percentage of cross-cluster individuals
total_indivs <- nrow(contemporary_HI)
cross_cluster <- 2/total_indivs

# sample size of mother-embryo sequencing
mother_embryo <- 49


# using rbinom()
# size = number of mother-embryo sequencing
# prob = probability of a cross-cluster individual
df <- rbinom(1000, mother_embryo, cross_cluster)

# the probability of finding AT LEAST 1 cross-cluster individual
find_prob <- sum(df>=1)/1000 # = 0.294
# the probability of NOT finding any cross-cluster individuals
not_find_prob <- 1-sum(df>=1)/1000 # = 0.706


#### plot the full results for the supplement ####
df_plotting <- as.data.frame(df)

ggplot(data=df_plotting) + 
  geom_bar(aes(x=df)) +
  xlab("# cross-cluster individuals found") +
  scale_x_continuous(breaks = c(0:5)) + 
  ylab("Count") +
  theme_bw() +
  theme(legend.position="none", axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"))
# export as 5" x 4" PDF

