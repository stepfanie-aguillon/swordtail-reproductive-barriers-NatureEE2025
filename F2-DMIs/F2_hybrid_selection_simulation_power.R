a=0.5
mm=a^2
mb=2*a*(1-a)
bb=(1-a)^2

simulation_results<-{}

N_indiv=163

for(x in 1:1000){

s=0.8
h=0.5

N=2000

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

quantile(simulation_results[,5])


thresh<-0.42

length(subset(simulation_results[,5],simulation_results[,5]<thresh))/1000

######results:
s_coeff<-c(0.3,0.4,0.5,0.6,0.7,0.8)
power<-c(0.08,0.24,0.55,0.83,0.99,1)