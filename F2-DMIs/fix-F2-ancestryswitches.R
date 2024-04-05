

#### getGenotypeSegments_fromGenotypes.R

options(stringAsFactors=FALSE)
library(data.table)

#args <- commandArgs(TRUE)
genotypeBase <- "XcorXbir-F2-fullgenotypes_selectshared_021924.txt"
#genotypeTransposedFile <- paste("./data/", genotypeBase, "_transposed", sep="")
genotypeTransposedFile <- "XcorXbir-F2-fullgenotypes_selectshared_021924.txt"


#if (!file.exists(genotypeTransposedFile)) {
#  cmd0=paste("perl /home/groups/schumer/shared_bin/Lab_shared_scripts/transpose_nameout.pl ",genotypeBase,sep="");
#  system(cmd0)
#}

genoAllInds <- fread(genotypeTransposedFile, header=T, sep="\t", fill=TRUE)
genoAllInds <- data.frame(genoAllInds)
genoAllInds$chrID <- as.character(data.frame(do.call('rbind', strsplit(as.character(genoAllInds$id), ":")))$X1)
genoAllInds$pos <- as.character(data.frame(do.call('rbind', strsplit(as.character(genoAllInds$id), ":")))$X2)
allchrs<-unique(genoAllInds$chr)

genoHeader <- colnames(genoAllInds)
allSeg <- data.frame()
for (i in 2:(length(genoHeader)-2)) {
  sampleID <- genoHeader[i]
  #print(sampleID)
  indGenoWG<-data.frame(chrPos=genoAllInds$id, chr=genoAllInds$chr, pos=as.numeric(genoAllInds$pos), geno=genoAllInds[,i])
  indGenoWG <- na.omit(indGenoWG)
  indRunsWG <- data.frame()
  
  for (k in 1:length(allchrs)){
    chrID <- as.character(allchrs[k])
    indChr<-subset(indGenoWG, chr==chrID)
    
    if (nrow(indChr)>0) {
      indGenoRuns<-{}
      first_geno=indChr$geno[1]
      start=indChr$pos[1]
      for(x in 1:nrow(indChr)){
        #print(indChr$pos[x])
        focal=indChr$geno[x]
        if(focal != first_geno){
          stop=indChr$pos[x-1]
          indGenoRuns<-rbind(indGenoRuns,cbind(start,stop,first_geno))
          first_geno=focal
          start=indChr$pos[x]
        }
        if (x==nrow(indChr)) {
          stop=indChr$pos[x]
          indGenoRuns<-rbind(indGenoRuns,cbind(start,stop,first_geno))
        }
      }
      if (!(is.null(indGenoRuns))) {
        indGenoRunsChr <- data.frame(chr=chrID, indGenoRuns)
        indRunsWG <- rbind(indRunsWG, indGenoRunsChr)
      }
    }
  }
  
  colnames(indRunsWG)[4] ="genotype"
  # write.table(indRunsWG, paste(sampleID, "_genotype_segments.txt", sep=""), row.names = F, quote = F, sep="\t")
  indRunsWG$indiv <- sampleID
  allSeg <- rbind(allSeg, indRunsWG)
  #hmzPar0 <- subset(indRunsWG,genotype==0)
  #hmzPar2 <- subset(indRunsWG,genotype==2)
  #het <- subset(indRunsWG,genotype==1)
  #write.table(hmzPar0, paste(sampleID, "_homozygous_parent0_segments.txt", sep=""), row.names = F, quote = F, sep="\t")
  #write.table(hmzPar2, paste(sampleID, "_homozygous_parent2_segments.txt", sep=""), row.names = F, quote = F, sep="\t")
  #write.table(het, paste(sampleID, "_heterozygous_segments.txt", sep=""), row.names = F, quote = F, sep="\t")
}
write.table(allSeg, paste(genotypeTransposedFile, "_allInds_genotype_segments.txt", sep=""), row.names = F, quote = F, sep="\t")








#### filterGenotypeSegments_getCrossovers.R

options(stringAsFactors=FALSE)
#args <- commandArgs(TRUE)
input <- "XcorXbir-F2-fullgenotypes_selectshared_021924.txt"
#thresh <- 250000
thresh <- 5000
#print(input)


#genoSegsWG <- read.table("genotypes_CHPL2021_corCluster.txt_transposed_allInds_genotype_segments.txt", header=T)
genoSegsWG <- read.table(paste(input, "_allInds_genotype_segments.txt", sep=""), header=T)

genoSegsWG$segLen <- genoSegsWG$stop - genoSegsWG$start
allinds <- unique(genoSegsWG$indiv)
genoSegsWG <- subset(genoSegsWG, segLen!=0)

genoSegsWG$filter <- NA
genoSegsWG$filter[which(genoSegsWG$segLen<=thresh)] <- "remove"
genoSegsWG$filter[which(genoSegsWG$segLen>thresh)] <- "pass"

allchrs <- unique(genoSegsWG$chr)

filteredDF <- data.frame()
for (chrID in allchrs) {
  #print(chrID)
  #print(which(allchrs==chrID))
  chrSegs <- subset(genoSegsWG, chr==chrID)
  for (indID in allinds) {
    #print(indID)
    indChrSegs <- subset(chrSegs, indiv==indID)
    filterIndChr <- data.frame()
    i=1
    while (i < nrow(indChrSegs)) {
      #print(i)
      if (indChrSegs$filter[i] == "pass" & indChrSegs$filter[i+1] == "pass") {
        filterIndChr <- rbind(filterIndChr, indChrSegs[i,])
        i <- i + 1
      } else if (indChrSegs$filter[i] == "pass" & indChrSegs$filter[i+1] == "remove") {
        if (i==nrow(indChrSegs)-1) {
          filterIndChr <- rbind(filterIndChr, indChrSegs[i,])
          i <- i + 1
        }
        j <- i + 2
        while (j <= nrow(indChrSegs)) {
          if (indChrSegs$filter[j] == "pass") {
            if (indChrSegs$geno[i] != indChrSegs$geno[j]) {
              filterIndChr <- rbind(filterIndChr, indChrSegs[i,], indChrSegs[j,])
              i <- j +1
              j <- nrow(indChrSegs)+1
            } else if (indChrSegs$geno[i] == indChrSegs$geno[j]) {
              newSeg <- data.frame(chr=chrID, start=indChrSegs$start[i], stop=indChrSegs$stop[j], genotype=indChrSegs$geno[i], indiv=indID, segLen=indChrSegs$stop[j]-indChrSegs$start[i], filter=NA)
              filterIndChr <- rbind(filterIndChr, newSeg)
              i <- j +1
              j <- nrow(indChrSegs)+1
            } else {
              print(paste("Check merging i =", i, " j =", j))
            }
          } else if (j==nrow(indChrSegs)) {
            filterIndChr <- rbind(filterIndChr, indChrSegs[i,])
            i <- j +1
            j <- nrow(indChrSegs)+1
          } else {
            j <- j + 1
          }
        }
      } else if (indChrSegs$filter[i] == "remove" & i!=1) {
        if (nrow(filterIndChr)>0) {
          preMerged <- filterIndChr[nrow(filterIndChr),]
          j <- i + 1
          while (j <= nrow(indChrSegs)) {
            if (indChrSegs$filter[j] == "pass") {
              if (preMerged$geno != indChrSegs$geno[j]) {
                filterIndChr <- rbind(filterIndChr, indChrSegs[j,])
                i <- j + 1
                j <- nrow(indChrSegs)+1
              } else if (preMerged$geno == indChrSegs$geno[j]) {
                newSeg <- data.frame(chr=chrID, start=preMerged$start, stop=indChrSegs$stop[j], genotype=preMerged$geno, indiv=indID, segLen=indChrSegs$stop[j]-preMerged$start, filter=NA)
                filterIndChr[nrow(filterIndChr),] <- newSeg
                i <- j +1
                j <- nrow(indChrSegs)+1
              } else {
                print(paste("Check merging i =", i, " j =", j))
                j <- j + 1
              }
            } else {
              j <- j + 1
              i = j
            }
          }
        } else {
          i <- i + 1
        }
      } else if (indChrSegs$filter[i] == "remove" & i==1) {
        i <- i + 1
      } else {
        print(paste("Check pass fail i =", i))
        i <- i + 1
      }
    }
    if (i == nrow(indChrSegs) & indChrSegs$filter[i]=="pass") {
      filterIndChr <- rbind(filterIndChr, indChrSegs[i,])
    }
    filterIndChr$distBwt <- filterIndChr$start-c(0,filterIndChr$stop[-nrow(filterIndChr)])
    filteredDF <- rbind(filteredDF, filterIndChr)
  }
}

write.table(filteredDF, paste(input, "_allInds_genotype_segments.filtered_5kb.txt", sep=""), row.names = F, quote=F, sep="\t")

crossSum <- data.frame()
crossDF <- data.frame()
for (chrID in allchrs) {
  #print(chrID)
  #print(which(allchrs==chrID))
  chrSegs <- subset(filteredDF, chr==chrID)
  for (indID in allinds) {
    #print(indID)
    indChrSegs <- subset(chrSegs, indiv==indID)
    if (nrow(indChrSegs)>0) {
      chrIndSum <- data.frame(chr=chrID, indiv=indID, nCrossovers=(nrow(indChrSegs)-1), medSegLen=median(indChrSegs$segLen), medDistBtw=median(indChrSegs$distBwt))
      crossSum <- rbind(crossSum, chrIndSum)
      if (nrow(indChrSegs)>1) {
        for (i in 1:(nrow(indChrSegs)-1)) {
          tmpSwitch <- data.frame(chrName=chrID, start=indChrSegs$stop[i], stop=indChrSegs$start[i+1], sampleID=indID, genoSwitch=paste(indChrSegs$genotype[i], indChrSegs$genotype[i+1], sep="-"))
          crossDF <- rbind(crossDF, tmpSwitch)
        }
      }
    }
  }
}

write.table(crossSum, paste(input, "_allInds_crossover_counts.chrSummary_5kb.txt", sep=""), row.names = F, quote=F, sep="\t")
write.table(crossDF, paste(input, "_allInds_crossover_intervals_5kb.txt", sep=""), row.names = F, quote=F, sep="\t")








###############
# starting to correct for this by identifying common start/stops at the segments

library(tidyverse)



### get start/stops for chromosome 13

filtered_segments <- read_tsv("./XcorXbir-F2-fullgenotypes_selectshared_021924.txt_allInds_genotype_segments.filtered_5kb.txt", col_names=TRUE) %>%
  filter(chr=="chr-13")

stop_group <- filtered_segments %>%
  group_by(stop) %>%
  summarize(count = n()) %>%
  filter(count > 1)
# 20538312
# 20868571


start_group <- filtered_segments %>%
  group_by(start) %>%
  summarize(count = n()) %>%
  filter(count > 1)
# 20751393
# 20875924

# filter line used for plotting: 
# filter(!between(position, 20538312, 20875924))



### get start/stops for chromosome 6

filtered_segments <- read_tsv("./XcorXbir-F2-fullgenotypes_selectshared_021924.txt_allInds_genotype_segments.filtered_5kb.txt", col_names=TRUE) %>%
  filter(chr=="chr-06")

stop_group <- filtered_segments %>%
  group_by(stop) %>%
  summarize(count = n()) %>%
  filter(count > 1)
# 7754852
# 7966895

# 23899215
# 24120326

start_group <- filtered_segments %>%
  group_by(start) %>%
  summarize(count = n()) %>%
  filter(count > 1)
# 7764777
# 8085477

# 23914586
# 24123182


# filter line used for plotting: 
# filter(!between(position, 7754852, 8085477)) %>%
# filter(!between(position, 23899215, 24123182))




### get start/stops for chromosome 7


filtered_segments <- read_tsv("./XcorXbir-F2-fullgenotypes_selectshared_021924.txt_allInds_genotype_segments.filtered_5kb.txt", col_names=TRUE) %>%
  filter(chr=="chr-07")

stop_group <- filtered_segments %>%
  group_by(stop) %>%
  summarize(count = n()) %>%
  filter(count > 1)
# 24792083
# 25721848


start_group <- filtered_segments %>%
  group_by(start) %>%
  summarize(count = n()) %>%
  filter(count > 1)
# 24846367
# 25737000


# filter line used for plotting: 
# filter(!between(position, 24792083, 25737000))





filtered_segments <- read_tsv("./XcorXbir-F2-fullgenotypes_selectshared_021924.txt_allInds_genotype_segments.filtered_5kb.txt", col_names=TRUE) %>%
  filter(chr=="chr-14")

stop_group <- filtered_segments %>%
  group_by(stop) %>%
  summarize(count = n()) %>%
  filter(count > 1)
# 18269172
# 18519754


start_group <- filtered_segments %>%
  group_by(start) %>%
  summarize(count = n()) %>%
  filter(count > 1)
# 18405284
# 18539743


# filter line used for plotting: 
# filter(!between(position, 18269172, 18539743))

