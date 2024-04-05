
###################
###LOAD PACKAGE
##################
#!install.packages("truncnorm")
library(truncnorm)

########################
#####SET PARAMETERS
#######################
n=17
#number of individuals

total_time=600
#time of the trial

mean_group1=235.05
#mean time spent with cue 1

mean_group2=97.69
#mean time spent with cue 2

sd_group1=127.94
#sd of time spent with cue 1

sd_group2=81.36
#sd of time spent with cue 2

#####################
####RUN SIMULATION
#####################
#initialize simulation
x=0
pvals<-{}
ratio<-{}
array1<-{}
array2<-{}
for (k in 1:1000){

while (x < n){
q1<-round(rtruncnorm(1,a=10,b=total_time,mean=mean_group1,sd=sd_group1))
q2<-round(rtruncnorm(1,a=10,b=total_time,mean=mean_group2,sd=sd_group2))
	if ((q1+q2)<total_time){
	x=x+1
	array1<-c(array1,q1)
	array2<-c(array2,q2)
	}#if a possible observation
}#until the simulation is accepted
pvals<-c(wilcox.test(array1,array2,paired=TRUE,exact=FALSE)$p.value,pvals)
ratio<-c(ratio,mean(array1/(array1+array2)))
x=0
array1<-{}
array2<-{}
}#for 1000 simulations

#power is given by:
length(subset(pvals,pvals<0.05))/length(pvals)

####################
###with weaker preferences
####################



###################
###LOAD PACKAGE
##################
#!install.packages("truncnorm")
library(truncnorm)

########################
#####SET PARAMETERS
#######################
n=18
#number of individuals

total_time=600
#time of the trial

mean_group1=235.05*0.8
#mean time spent with cue 1

mean_group2=97.69*1.2
#mean time spent with cue 2

sd_group1=127.94
#sd of time spent with cue 1

sd_group2=81.36
#sd of time spent with cue 2

#####################
####RUN SIMULATION
#####################
#initialize simulation
x=0
pvals<-{}
ratio<-{}
array1<-{}
array2<-{}
for (k in 1:1000){

while (x < n){
q1<-round(rtruncnorm(1,a=10,b=total_time,mean=mean_group1,sd=sd_group1))
q2<-round(rtruncnorm(1,a=10,b=total_time,mean=mean_group2,sd=sd_group2))
	if ((q1+q2)<total_time){
	x=x+1
	array1<-c(array1,q1)
	array2<-c(array2,q2)
	}#if a possible observation
}#until the simulation is accepted
pvals<-c(wilcox.test(array1,array2,paired=TRUE,exact=FALSE)$p.value,pvals)
ratio<-c(ratio,mean(array1/(array1+array2)))
x=0
array1<-{}
array2<-{}
}#for 1000 simulations

#power is given by:
length(subset(pvals,pvals<0.05))/length(pvals)
