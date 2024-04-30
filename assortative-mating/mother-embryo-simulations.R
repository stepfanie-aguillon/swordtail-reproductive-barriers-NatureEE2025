# This script contains the R code to process mother/embryo data and run assortative mating simulations
#
# Authors: Aguillon SM, et al.
# Year: 2024
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: TBD
# bioRxiv DOI: https://doi.org/10.1101/2024.04.16.589374
#
# Edited date: 5 Oct 2023
#
# Please cite the paper if you use these scripts
#


# load packages
library(Rmisc)
library(tidyverse) #need to load AFTER Rmisc, so some functions work correctly
library(cowplot)

# load data
contemporary_HI <- read_tsv("./data/CHPL-contemporary-ancestry.txt", col_names=TRUE)
mother_embryo <- read_tsv("./data/mother-embyro-ancestry.txt", col_names=TRUE)

# plot colors
corcol <- "#359D6F"
bircol <- "#3A6BC6"



##### DATA PROCESSING #####
# need to ID mothers/embryos and match appropriately

# ID embryos
mother_embryo <- mother_embryo %>%
  mutate(age = ifelse(str_detect(sample_ID, "-E"), "embryo", "adult")) %>%
  mutate(matching_ID = sample_ID)
# update sample naming for matching between mother and embryo
#mothers with IDs like: CHPL-3-VI-21-222 --> CHPL-3-VI-21-F222 (n=8) 
mother_embryo$matching_ID <- str_replace(mother_embryo$matching_ID, "CHPL-3-VI-21-", "CHPL-3-VI-21-F")
#embryos with IDs like: CHPL-VI-21-F222-E2 --> CHPL-3-VI-21-F222-E2 (n=20)
mother_embryo$matching_ID <- str_replace(mother_embryo$matching_ID, "CHPL-VI-21-F", "CHPL-3-VI-21-F")
#samples with IDs like: CHPL-21-III-22-F-74 --> CHPL-21-III-22-F74 (n=34)
mother_embryo$matching_ID <- str_replace(mother_embryo$matching_ID, "CHPL-21-III-22-F-", "CHPL-21-III-22-F")
#embryos with IDs like: CHPL-21-III-22-38F --> CHPL-21-III-22-F38 (n=22)
mother_embryo$matching_ID <- str_replace(mother_embryo$matching_ID, "CHPL-21-III-22-38F", "CHPL-21-III-22-F38")
mother_embryo$matching_ID <- str_replace(mother_embryo$matching_ID, "CHPL-21-III-22-43F", "CHPL-21-III-22-F43")
mother_embryo$matching_ID <- str_replace(mother_embryo$matching_ID, "CHPL-21-III-22-49F", "CHPL-21-III-22-F49")
mother_embryo$matching_ID <- str_replace(mother_embryo$matching_ID, "CHPL-21-III-22-56F", "CHPL-21-III-22-F56")
mother_embryo$matching_ID <- str_replace(mother_embryo$matching_ID, "CHPL-21-III-22-57F", "CHPL-21-III-22-F57")
mother_embryo$matching_ID <- str_replace(mother_embryo$matching_ID, "CHPL-21-III-22-58F", "CHPL-21-III-22-F58")
mother_embryo$matching_ID <- str_replace(mother_embryo$matching_ID, "CHPL-21-III-22-60F", "CHPL-21-III-22-F60")
mother_embryo$matching_ID <- str_replace(mother_embryo$matching_ID, "CHPL-21-III-22-66F", "CHPL-21-III-22-F66")
mother_embryo$matching_ID <- str_replace(mother_embryo$matching_ID, "CHPL-21-III-22-67F", "CHPL-21-III-22-F67")
mother_embryo$matching_ID <- str_replace(mother_embryo$matching_ID, "CHPL-21-III-22-70F", "CHPL-21-III-22-F70")
mother_embryo$matching_ID <- str_replace(mother_embryo$matching_ID, "CHPL-21-III-22-73F", "CHPL-21-III-22-F73")
#remove REDO and .fastq from endings
mother_embryo$matching_ID <- str_remove(mother_embryo$matching_ID, ".R1.fastq")
mother_embryo$matching_ID <- str_remove(mother_embryo$matching_ID, "_REDO")


# separate embryos from adults and edit dataset to make merging easier
# final df has one column for mother's ID and one column for embryo number
embryos_df <- mother_embryo %>%
  filter(age=="embryo") %>%
  select(matching_ID, cor_count, bir_count, hybrid_index, heterozygosity, read_count) %>%
  rename(embryo_ID = matching_ID)
# to get splitting to work correctly, need to slightly modify names (will edit back after)
embryos_df$embryo_ID <- str_replace(embryos_df$embryo_ID, "CHPL-V-22", "CHPL-1-V-22")
embryos_df <- embryos_df %>%
  separate(embryo_ID, c("A","B","C","D","E","F"), remove=FALSE) %>%
  unite(ID, A, B, C, D, E, sep="-") %>%
  dplyr::rename(embryo_num = F, embryo_corcount = cor_count, embryo_bircount = bir_count, embryo_hybrid_index = hybrid_index, embryo_heterozygosity = heterozygosity, embryo_read_count = read_count) %>%
  select(ID:embryo_read_count)
# fix names
embryos_df$ID <- str_replace(embryos_df$ID, "CHPL-1-V-22", "CHPL-V-22")


# separate adults from embryos
adults_df <- mother_embryo %>%
  filter(age=="adult") %>%
  select(matching_ID, cor_count, bir_count, hybrid_index, heterozygosity, read_count) %>%
  rename(ID = matching_ID) %>%
  mutate(ancestry_group = ifelse(hybrid_index > 0.6, "cortezi", ifelse(hybrid_index < 0.1, "birchmanni", "cross-cluster")))


# join separated mother and embryo datasets
#inner_join() to make sure only individuals that match are included
matched_mother_embryo_df <- inner_join(embryos_df, adults_df, by="ID")

# calculate differences between mother and embryo ancestries
matched_mother_embryo_df <- matched_mother_embryo_df %>%
  mutate(diff = hybrid_index - embryo_hybrid_index) %>%
  drop_na() %>%
  arrange(hybrid_index) %>%
  mutate(plot_order=1:nrow(matched_mother_embryo_df))







##### PREPPING SIMULATIONS #####

# full dataset of observed mother-embryo sequencing
# most females have multiple embryos (depending on sequencing success)
# plot to check for any large amount of within-individual variation
ggplot(data=matched_mother_embryo_df) + 
  geom_vline(aes(xintercept=0), color="dark gray", linewidth=1) + 
  geom_point(aes(x=diff, y=plot_order, color=ancestry_group), size=2) +
  #facet_grid(~ancestry_group, scales="free") +
  scale_color_manual(values=c(bircol, corcol), name="Maternal ancestry", labels=c("birchmanni cluster","cortezi-like cluster")) +
  xlab("Maternal/Embryo ancestry difference") +
  xlim(c(-0.5,0.5)) +
  theme_bw() +
  theme(legend.position=c(0.2,0.8),legend.margin=margin(-2),axis.line.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank()) +
  theme(axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"), strip.text.x=element_text(size=12,face="bold"))
# export as 5" x 4" PDF
# in Illustrator: italicized species names, fixed look of legend


# largely following script available here: https://github.com/Schumerlab/Xbir_xcor_hybridzone/blob/master/assortative_mating_analysis.R

# subset dataset to keep only one mother-offspring pair per individual
matched_subset <- matched_mother_embryo_df %>%
  group_by(ID) %>%
  sample_n(size=1)

# observed maternal-offspring difference in ancestry_group overall
observed_diff <- mean(abs(matched_subset$embryo_hybrid_index - matched_subset$hybrid_index))
observed_diff_CI <- CI(abs(matched_subset$embryo_hybrid_index - matched_subset$hybrid_index),ci=0.95)



### PARAMETERS ###

###CLUSTER THRESHOLD = maximum difference in ancestry within a cluster

# using a random normal distribution to estimate this based on mean and sd of the observed data
threshold_birchmanni <- contemporary_HI %>%
  filter(hybrid_index<0.1) %>%
  summarize(mean = mean(hybrid_index), sd = sd(hybrid_index))

threshold_cortezi <- contemporary_HI %>%
  filter(hybrid_index>0.6) %>%
  summarize(mean = mean(hybrid_index), sd = sd(hybrid_index))

rnorm_cortezi <- rnorm(10000,mean=threshold_cortezi$mean,sd=threshold_cortezi$sd)
rnorm_birchmanni <- rnorm(10000,mean=threshold_birchmanni$mean,sd=threshold_birchmanni$sd)

diff_cortezi <- ifelse(max(rnorm_cortezi)>1,1,max(rnorm_cortezi)) - min(rnorm_cortezi)
diff_birchmanni <- max(rnorm_birchmanni) - ifelse(min(rnorm_birchmanni)<0,0,min(rnorm_birchmanni))

cluster_threshold <- round(max(diff_cortezi,diff_birchmanni),4)


###ANCESTRY SD = variance in ancestry between siblings

# calculating from (full) observed data
ancestry_sd <- matched_mother_embryo_df %>%
  group_by(ID) %>%
  summarize(offspring_ancestry_sd = sd(embryo_hybrid_index))

ancestry_sd <- round(mean(na.omit(ancestry_sd$offspring_ancestry_sd)),4)


###ADDITIONAL PARAMETERS
number_of_simulations <- 500
sample_size <- as.numeric(nrow(matched_subset)) #number of mother-offspring pairs to simulate




##### RUNNING SIMULATIONS #####

# vary assortative mating value over the interval of 0 to 1
mating_prop = seq(0,1,by=0.01)
# or set to a specific value (e.g., to simulate COMPLETE or RANDOM assortative mating)
#mating_prop = 1 #COMPLETE
#mating_prop = 0 #RANDOM

# make output matrix for storing results
output <- matrix(NA,nrow=length(mating_prop),ncol=5)


# simulation: loops over the mating proportion interval with 500 simulated matings for each value
for (P in 1:length(mating_prop)) {
  mean <- {}
  
  for (j in 1:number_of_simulations){
    mates<-c(adults_df$hybrid_index)
    offspring_index<-{}
    maternal_index<-{}
    for (x in 1:length(matched_subset$hybrid_index)){
      maternal<-matched_subset$hybrid_index[x]
      mate<-sample(mates,1)
      hold=0
      while(hold==0){
        if(abs(maternal-mate)<cluster_threshold){
          hold=1
        }
        if(abs(maternal-mate)>cluster_threshold){
          pref=rbinom(1,1,1-mating_prop[P])
          if(pref == 1){
            hold=1
          } else{
            mate<-sample(mates,1)
          }
        }
        
      }
      
      off<-rnorm(1,mean=(maternal+mate)/2,sd=ancestry_sd)
      if(off<0){off=0}#don't allow hybrid index <1
      if(off>1){off=1}#don't allow hybrid index >1
      offspring_index<-c(offspring_index,off)
      maternal_index<-c(maternal_index,maternal)
      
    }
    new<-cbind(offspring_index,maternal_index)
    new<-new[sample(nrow(new),sample_size),]
    mean<-c(mean,mean(abs(new[,1]-new[,2])))
    
  }
  
  mean_assort <- mean
  mean_assort_CI <- CI(mean_assort,ci=0.95)
  
  ### summarize output
  # mating prop value in simulation
  output[P,1] <- mating_prop[P]
  
  # confidence intervals for difference between mother and offspring value
  output[P,2] <- mean_assort_CI[1] #lower 95%
  output[P,3] <- mean_assort_CI[2] #mean
  output[P,4] <- mean_assort_CI[3] #upper 95%
  
  # proportion of simulations with stronger assortative mating than observed data
  output[P,5] <- length(subset(mean_assort,mean_assort<observed_diff))/length(mean_assort)
  
}

# rerun simulations with specific parameters and save as follows:
#write.table(new, "./100-assortative-mating-sim.txt", quote=FALSE, sep="\t", row.names=FALSE)
#write.table(new, "./random-mating-sim.txt", quote=FALSE, sep="\t", row.names=FALSE)

#### PLOT RESULTS ####

output <- as.data.frame(output) %>%
  rename(assort_mating_value = V1,
         sim_lowerCI = V2,
         sim_mean = V3,
         sim_upperCI = V4,
         prop_stronger_sims = V5)

output <- output %>%
  mutate(sim_diff = sim_mean - observed_diff_CI[2])

#### model selection plot for supplement #### 
ggplot() +
  geom_hline(aes(yintercept=0),color="dark gray",linewidth=1) +
  geom_smooth(data=output,aes(x=assort_mating_value,y=sim_diff)) +
  geom_point(data=output,aes(x=assort_mating_value,y=sim_diff),shape=21) +
  xlab("Simulated assortative mating value") +
  ylab("Simulated mean - Observed mean") +
  theme_bw() +
  theme(legend.position="none", axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))
# export as 5" x 4" PDF

#100% assortative mating minimizes the difference between the simulated and observed datasets





##### SIMULATION RESULTS #####

matched_subset <- matched_subset %>%
  ungroup() %>%
  arrange(hybrid_index) %>%
  mutate(plot_order = 1:nrow(matched_subset), dataset="observed") 

### NOTE!!!: must first run simulations for specific assortative mating values and save results (see lines 177 and 244)

# load and wrangle simulation results
assort_sim <- read_tsv("100-assortative-mating-sim.txt", col_names=TRUE)
random_sim <- read_tsv("random-mating-sim.txt", col_names=TRUE)

assort_sim <- assort_sim %>%
  mutate(diff = maternal_index - offspring_index) %>%
  mutate(ancestry_group = ifelse(maternal_index > 0.5,"cortezi","birchmanni")) %>%
  arrange(maternal_index) %>%
  mutate(plot_order = 1:nrow(assort_sim), dataset = "simulation") %>%
  rename(hybrid_index_mother = maternal_index, hybrid_index_embryo = offspring_index)

random_sim <- random_sim %>%
  mutate(diff = maternal_index - offspring_index) %>%
  mutate(ancestry_group = ifelse(maternal_index > 0.5,"cortezi","birchmanni")) %>%
  arrange(maternal_index) %>%
  mutate(plot_order = 1:nrow(random_sim))

## estimation of difference between mother and embryo with cross-cluster mating
random_sim_filter <- random_sim %>%
  mutate(abs_diff = abs(diff)) %>%
  filter(abs_diff > 0.05) %>%
  summarize(n=n(),
            avg_abs_diff = mean(abs_diff),
            sd_abs_diff = sd(abs_diff))


#### plotting observed data with 100% assortative mating simulations (Figure 2A) ####
panel1 <- ggplot() + 
  geom_vline(aes(xintercept=0), color="dark gray", size=1) + 
  geom_point(data=assort_sim, aes(x=diff, y=plot_order, color=ancestry_group), size=2.5, shape=21) + 
  geom_point(data=matched_subset, aes(x=diff, y=plot_order, color=ancestry_group), size=2.5) +
  scale_color_manual(values=c(bircol,corcol),name="Maternal ancestry",labels=c("birchmanni cluster","cortezi-like cluster")) +
  xlab("") +
  xlim(c(-0.4,0.4)) +
  theme_bw() +
  theme(legend.position=c(0.22,0.8),legend.margin=margin(-2),axis.line.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank()) +
  theme(axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=8,color="black"), legend.title=element_text(size=8), legend.text=element_text(size=8))
# March 2024 note: this ended up being combined with behavior results

#### plotting random mating simulations ####
panel2 <- ggplot() + 
  geom_vline(aes(xintercept=0), color="dark gray", size=1) + 
  geom_point(data=random_sim, aes(x=diff, y=plot_order, color=ancestry_group), size=2.5, shape=2) + 
  scale_color_manual(values=c(bircol,corcol)) +
  xlab("Maternal/embryo ancestry difference") +
  xlim(c(-0.5,0.5)) +
  theme_bw() +
  theme(legend.position="none",legend.margin=margin(-2),axis.line.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank()) +
  theme(axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"))
# March 2024 note: this ended up being placed in the supplement


#### combine plots for in-text figure ####
plot_grid(panel1, panel2, align="v", nrow=2)
# export as 3.5" x 5" PDF
# in Illustrator: added plot labels, italicized species names, updated legend, fixed plot spacing
# March 2024 note: ended up combining `panel1` with behavior trial results and moving `panel2` to the supplement (exporting `panel2` in the same way as the full dataset fig above)
