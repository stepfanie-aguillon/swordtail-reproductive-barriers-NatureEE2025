getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

a=0.5
mm=a^2
mb=2*a*(1-a)
bb=(1-a)^2

simulation_results<-{}

for(x in 1:500000){

s=runif(1,0,1)
h=runif(1,0,1)

N=10000

fmm=mm*N
fmb=mb*N*(1-h*s)
fbb=bb*N*(1-s)

N_total=fmm+fmb+fbb

ratio_bb=fbb/N_total
ratio_mb=fmb/N_total
ratio_mm=fmm/N_total

pop<-c(rep(2,fmm),rep(1,fmb),rep(0,fbb))

N_indiv=163

sim_pop=sample(pop,N_indiv)

obb<-length(subset(sim_pop,sim_pop==0))
obm<-length(subset(sim_pop,sim_pop==1))

a<- (obb + 0.5*obm)/N_indiv

simulation_results<-rbind(simulation_results,cbind(s,h,obb,obm,a))

}


modes<-{}
for (k in 1:1000){

full<-c(rep(0,28),rep(1,587),rep(2,317))
sim<-sample(full,163)
sim_bb<-length(subset(sim,sim==0))
sim_bm<-length(subset(sim,sim==1))
sim_mm<-length(subset(sim,sim==2))
anc_sample<-1-mean(sim)/2

accepted<-subset(simulation_results,simulation_results[,5] >(anc_sample*0.95) & simulation_results[,5] <(anc_sample*1.05) & simulation_results[,3] <= round(sim_bb*1.05) & simulation_results[,3] >= round(sim_bb*0.95) & simulation_results[,4] <= round(sim_bm*1.05) & simulation_results[,4]>= round(sim_bm*0.95))
#all +/- binomial standard error

modes<-rbind(modes,cbind(getmode(round(accepted[,1],3)),getmode(round(accepted[,2],3)),quantile(accepted[,1],0.025),quantile(accepted[,1],0.975)))

}

colnames(modes)<-c("map_s","map_h","lower_bound_s","upper_bound_s")
subset(modes,modes[,3]<=0.69)

write.table(accepted,file="~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Mini_mito_project/Data/simulate_ndufa13_subsample163_range.txt",sep="\t",row.names=FALSE,quote=FALSE)