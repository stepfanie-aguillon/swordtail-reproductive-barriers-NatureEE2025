# This script contains the R code to plot F1 and F2 ancestry tracts
#
# Authors: Aguillon SM, et al.
# Year: 2025
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: Nature Ecology and Evolution
# DOI: https://doi.org/10.1038/s41559-025-02669-9
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
## genotypes for F2s at particular markers
F1_df <- read_tsv("./data/for-ancestry-tracts/XcorXbir-F1-repeats-genotypes.txt_transposed", col_names=TRUE)
PF2_df <- read_tsv("./data/for-ancestry-tracts/XcorXbir-PF2-genotypes-110623.txt_transposed", col_names=TRUE)
PF2_df2 <- read_tsv("./data/for-ancestry-tracts/XcorXbir-PF2-I-24-S255-genotypes-021924.txt_transposed", col_names=TRUE)


### plot example ancestry for F1s at two chromosomes
# chromosome 1
F1_df_filter <- F1_df %>%
  filter(chr=="chr-01") %>%
  select(chr,pos,"XcorXbir-PF1-B8-III-22-1-L-GGG-F.R1.fastq") %>%
  drop_na()
plot(main="XcorXbir-F1-B8-III-22-1-L-GGG-F", x=F1_df_filter$pos/1e6, y=F1_df_filter$`XcorXbir-PF1-B8-III-22-1-L-GGG-F.R1.fastq`, type="l", lwd=2, col="black", ylim=c(0,2), yaxt="n", xlab="Chromosome 1 Position (Mb)", ylab="# cortezi alleles")
axis(2, at=c(0,1,2))

# chromosome 2
F1_df_filter <- F1_df %>%
  filter(chr=="chr-02") %>%
  select(chr,pos,"XcorXbir-PF1-B8-III-22-1-L-GGG-F.R1.fastq") %>%
  drop_na()
plot(main="XcorXbir-F1-B8-III-22-1-L-GGG-F", x=F1_df_filter$pos/1e6, y=F1_df_filter$`XcorXbir-PF1-B8-III-22-1-L-GGG-F.R1.fastq`, type="l", lwd=2, col="black", ylim=c(0,2), yaxt="n", xlab="Chromosome 2 Position (Mb)", ylab="# cortezi alleles")
axis(2, at=c(0,1,2))

# export individual plots as 3.25" x 3.25" PDFs
# in Illustrator: combined into a single multi-panel, italicized species names, fixed plot spacing




### assess ancestry for F2 vs BC
# need to ensure that F2s are ACTUAL F2s (due to the possibility of retained sperm)

# make empty results df
# need data to be a data.frame for this to work !!
results <- {}
#data <- as.data.frame(PF2_df)
data <- as.data.frame(PF2_df2)

# calculation to assess 
for(x in 1:(length(data[1,])-2)){
  
  prop_bir<-length(subset(data[,(x+2)],data[,(x+2)]==0))/length(subset(data[,(x+2)],is.na(data[,(x+2)])==FALSE))
  
  prop_het<-length(subset(data[,(x+2)],data[,(x+2)]==1))/length(subset(data[,(x+2)],is.na(data[,(x+2)])==FALSE))
  
  prop_cor<-length( subset(data[,(x+2)],data[,(x+2)]==2))/length(subset(data[,(x+2)],is.na(data[,(x+2)])==FALSE))
  
  name<-colnames(data[x+2])
  
  results<-rbind(results,cbind(name,prop_bir,prop_het,prop_cor))
}
# save results
#write.table(results, "F2-vs-BC_comparison_110923.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(results, "F2-vs-BC_comparison_021924.txt", quote=FALSE, sep="\t", row.names=FALSE)


#individuals that are likely to be BC to cortezi
likely_BC <- subset(results, results[,2] < 0.05)





### plot example ancestry for F2s at two chromosomes
# chromosome 1
F2_df_filter <- PF2_df %>%
  filter(chr=="chr-01") %>%
  select(chr,pos,"XcorXbir-PF2-IX-23-S221-L-RRR-U-60.R1.fastq") %>%
  drop_na()
plot(main="XcorXbir-F2-IX-23-S221-L-RRR-U-60", x=F2_df_filter$pos/1e6, y=F2_df_filter$`XcorXbir-PF2-IX-23-S221-L-RRR-U-60.R1.fastq`, type="l", lwd=2, col="black", ylim=c(0,2), yaxt="n", xlab="Chromosome 1 Position (Mb)", ylab="# cortezi alleles")
axis(2, at=c(0,1,2))

# chromosome 2
F2_df_filter <- PF2_df %>%
  filter(chr=="chr-01") %>%
  select(chr,pos,"XcorXbir-PF2-IX-23-S221-L-YOO-U-19.R1.fastq") %>%
  drop_na()
plot(main="XcorXbir-F2-IX-23-S221-L-YOO-U-19", x=F2_df_filter$pos/1e6, y=F2_df_filter$`XcorXbir-PF2-IX-23-S221-L-YOO-U-19.R1.fastq`, type="l", lwd=2, col="black", ylim=c(0,2), yaxt="n", xlab="Chromosome 1 Position (Mb)", ylab="# cortezi alleles")
axis(2, at=c(0,1,2))

# export individual plots as 3.25" x 3.25" PDFs
# in Illustrator: combined into a single multi-panel, italicized species names, fixed plot spacing
