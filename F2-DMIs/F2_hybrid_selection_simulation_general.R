a=0.5
mm=a^2
mb=2*a*(1-a)
bb=(1-a)^2

simulation_results<-{}

#####chr7 parameters
chr="ch7"
focal_anc<-0.595092
focal_anc<- 1- focal_anc
obs_bb<-25
obs_cc<-56
obs_bc<-82

N_indiv=obs_bb+obs_cc+obs_bc

######chr14 parameters
chr="ch14"
focal_anc<-0.5828221
focal_anc<- 1- focal_anc
obs_bb<-35
obs_cc<-62
obs_bc<-66

N_indiv=obs_bb+obs_cc+obs_bc

#####chr13 parameters
chr="ch13"
focal_anc<-0.653
focal_anc<- 1- focal_anc
obs_bb<-0
obs_cc<-50
obs_bc<-113

N_indiv=obs_bb+obs_cc+obs_bc

#####chr6 parameters
chr="ch6"
focal_anc<-0.567
focal_anc<- 1- focal_anc
obs_bb<-23
obs_cc<-45
obs_bc<-95

N_indiv=obs_bb+obs_cc+obs_bc

#####set tolerance parameters for rejection sampling

tol=0.05
tol_up=1+tol
tol_down=1-tol

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


sim_pop=sample(pop,N_indiv)

obb<-length(subset(sim_pop,sim_pop==0))
obm<-length(subset(sim_pop,sim_pop==1))

a<- (obb + 0.5*obm)/N_indiv

simulation_results<-rbind(simulation_results,cbind(s,h,obb,obm,a))

}



accepted<-subset(simulation_results,simulation_results[,5] <(focal_anc*tol_up) & simulation_results[,5] > (focal_anc*0.95) & simulation_results[,3]>=obs_bb*tol_down & simulation_results[,3]<=obs_bb*tol_up & simulation_results[,4]>=obs_bc*tol_down & simulation_results[,4]<=obs_bc*tol_up)
length(accepted[,1])
#all +/- binomial standard error

write.table(accepted,file=paste("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Xcortezi_Xbirchmanni_reproductive_barriers/Data/ABC_simulations_accepted_F2s_Feb2024_",chr,".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)