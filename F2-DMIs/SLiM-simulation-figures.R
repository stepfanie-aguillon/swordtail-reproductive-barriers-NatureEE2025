# This script contains the R code to visualize SLiM simulation results
#
# Authors: Aguillon SM, et al.
# Year: 2024
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: TBD
# bioRxiv DOI: https://doi.org/10.1101/2024.04.16.589374
#
# Edited date: 2 Dec 2024
#
# Please cite the paper if you use these scripts
#

# load packages
library(tidyverse)
library(cowplot)

# set working directory
setwd("/Users/stepfanie/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Xcortezi_Xbirchmanni_reproductive_barriers/Data/mitochondrial_divergence_simulated")

# load data
sim_xmal_xcor <- read_csv("all_xmal_xcor_simulated_results_100sim_burnin.txt", col_names=FALSE)
sim_xmal_xmoz <- read_csv("all_xmal_xmoz_simulated_results_100sim_burnin.txt", col_names=FALSE)
sim_xmal_xbir <- read_csv("all_xmal_xbir_simulated_results_100sim_burnin.txt", col_names=FALSE)

# empirical averages
avg_xmal_xmal=0.008534004
avg_xcor_xcor=0.004164193
avg_xmal_xcor=0.008862201
avg_xmal_xmoz=0.035
avg_xmal_xbir=0.054

# plot colors
corcol <- "#359D6F"
bircol <- "#3A6BC6"
malcol <- "#AD0000"


### make density plots from SLiM results (Figure 3C) ###

# malinche versus cortezi simulations
panel1 <- ggplot(data=sim_xmal_xcor, aes(x=X1/16639)) + 
  geom_density(fill="dark gray", color="dark gray") +
  geom_vline(xintercept=avg_xmal_xmal, linetype="dashed", linewidth=1) +
  geom_vline(xintercept=avg_xcor_xcor, color=corcol, linewidth=1) + 
  geom_vline(xintercept=avg_xmal_xcor, color=malcol, linewidth=1) +
  xlim(0, 0.04) +
  theme_bw() +
  xlab("Simulated divergence: X. malinche - X. cortezi") + 
  ylab("Density") + 
  theme(legend.position="none", axis.title=element_text(face="bold",size=10), axis.text=element_text(size=8,color="black"),strip.text.x=element_text(size=8,face="bold"))

# malinche versus montezumae simulations
panel2 <- ggplot(data=sim_xmal_xmoz, aes(x=X1/16639)) + 
  geom_density(fill="dark gray", color="dark gray") +
  geom_vline(xintercept=avg_xmal_xmoz, linewidth=1) +
  xlim(0, 0.04) +
  theme_bw() +
  xlab("Simulated divergence: X. malinche - X. montezumae") + 
  ylab("Density") + 
  theme(legend.position="none", axis.title=element_text(face="bold",size=10), axis.text=element_text(size=8,color="black"),strip.text.x=element_text(size=8,face="bold"))

#combine plots
plot_grid(panel1, panel2, align="h", nrow=1)
# export as 5" tall x 7" wide
# in Illustrator: combine with other panels, resize and fix size of fonts, add in legend and panel labels




### make density plot from SLiM results for supplement ###

# malinche versus birchmanni simulations
ggplot(data=sim_xmal_xbir, aes(x=X1/16639)) + 
  geom_density(fill="dark gray", color="dark gray") +
  geom_vline(xintercept=avg_xmal_xbir, linewidth=1) +
  xlim(0, 0.06) +
  theme_bw() +
  xlab("Simulated divergence: X. malinche - X. birchmanni") + 
  ylab("Density") + 
  theme(legend.position="none", axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"),strip.text.x=element_text(size=12,face="bold"))
# export as 4" tall x 5" wide
# in Illustrator: italicize species names







#### revisions: divergence time estimate #### 

# load data
mt_div <- read_delim("./data/ABC_mito_tdiv_simulations_accepted_November2024.txt", delim="\t", col_names=TRUE)


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


## supplemental figure of results
ggplot() + 
  geom_density(data=mt_div, aes(x=Generations, y=after_stat(density*nrow(mt_div))*1000), fill="darkgray", color=NA) +
  geom_vline(xintercept=posterior.mode(mt_div$Generations), color="black", linetype="dashed") +
  ylab("Density") +
  xlab("Estimated divergence time (generations)") +
  theme_bw() +
  theme(legend.position="none", axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=10,color="black"))
# export as 4" high x 5" wide

posterior.mode(mt_div$Generations) #202320.2
round(quantile(mt_div$Generations,c(0.025,0.975)),3)
#2.5%    97.5% 
#165287.9 231917.5 
