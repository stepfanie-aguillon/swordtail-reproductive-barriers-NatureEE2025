file="XcorXbir-F2-fullgenotypes_selectshared_021924.txt_thinned_physical_dist.txt"
data<-read.csv(file=file,sep="\t",head=TRUE)


for(k in 1:500){

index<-read.csv(file="XcorXbir-F2-fullgenotypes_selectshared_021924.txt_hybrid_index",sep="\t",head=TRUE)
focal<- rbinom(length(index$hybrid_index),1,index$hybrid_index) + rbinom(length(index$hybrid_index),1,index$hybrid_index) #simulate focal locus

whole_genome<-{}

y=2
focal1<-focal

results<-{}


for(i in y:length(data[1,])){

focal2<-data[,i]

s1<-as.numeric(focal1[1:length(data[,1])])
s2<-as.numeric(focal2[1:length(data[,1])])

combined<-na.omit(cbind(s1,s2))
total<-length(combined[,1])
BB<-length(subset(combined[,1],combined[,1]==0))/total
MB<-length(subset(combined[,1],combined[,1]==1))/total
MM<-length(subset(combined[,1],combined[,1]==2))/total

bb<-length(subset(combined[,1],combined[,2]==0))/total
mb<-length(subset(combined[,1],combined[,2]==1))/total
mm<-length(subset(combined[,1],combined[,2]==2))/total

exp1<-MM*mm*total
exp2<-MM*mb*total
exp3<-MM*bb*total
exp4<-MB*mm*total
exp5<-MB*mb*total
exp6<-MB*bb*total
exp7<-BB*mm*total
exp8<-BB*mb*total
exp9<-BB*bb*total

obs1<-length(subset(combined[,1],(combined[,1]==2 & combined[,2]==2) ))
obs2<-length(subset(combined[,1],(combined[,1]==2 & combined[,2]==1) ))
obs3<-length(subset(combined[,1],(combined[,1]==2 & combined[,2]==0) ))
obs4<-length(subset(combined[,1],(combined[,1]==1 & combined[,2]==2) ))
obs5<-length(subset(combined[,1],(combined[,1]==1 & combined[,2]==1) ))
obs6<-length(subset(combined[,1],(combined[,1]==1 & combined[,2]==0) ))
obs7<-length(subset(combined[,1],(combined[,1]==0 & combined[,2]==2) ))
obs8<-length(subset(combined[,1],(combined[,1]==0 & combined[,2]==1) ))
obs9<-length(subset(combined[,1],(combined[,1]==0 & combined[,2]==0) ))

chival<- (((obs1-exp1)^2)/exp1) + (((obs2-exp2)^2)/exp2) + (((obs3-exp3)^2)/exp3) + (((obs4-exp4)^2)/exp4) + (((obs5-exp5)^2)/exp5) + (((obs6-exp6)^2)/exp6) + (((obs7-exp7)^2)/exp7) + (((obs8-exp8)^2)/exp8) + (((obs9-exp9)^2)/exp9)

results<-rbind(results,cbind("null_peak",as.character(colnames(data[i])),mean(combined[,1])/2,mean(combined[,2])/2,exp1,exp2,exp3,exp4,exp5,exp6,exp7,exp8,exp9,obs1,obs2,obs3,obs4,obs5,obs6,obs7,obs8,obs9,total,chival))

}

whole_genome<-rbind(whole_genome,results)

write.table(whole_genome,file=paste("results_null_chisq_chr7_chr14_genome_scan_",file,"_",k,sep=""),sep="\t",row.names=FALSE,col.names=c("chr1_site1","chr2_site2","a1","a2","exp1","exp2","exp3","exp4","exp5","exp6","exp7","exp8","exp9","MMmm","MMmb","MMbb","MBmm","MBmb","MBbb","BBmm","BBmb","BBbb","total","chisqstat"),quote=FALSE)


}
