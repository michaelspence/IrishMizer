setwd("C:/Users/KA03/CEFAS/Michael Spence (Cefas) - Inference on a budget/mizer/IrishSea")
library(mizer)

#source("MixedFish/multiple_gears.R")

dat<-read.csv("Data/IS_Params_7.csv")
dat<-rbind(dat,dat[4,])
dat$species <- as.character(dat$species)
dat$species[4]<-"Nephrops_w"
dat$species[8]<-"Nephrops_e"
dat$W_a[4] <- dat$W_a[8] <- dat$W_a[4] * 1e3 ### adjustment for units
dat$w_inf<-dat$W_a*(dat$Linf^dat$W_b)
dat$w_mat<-dat$W_a*(dat$Lmat^dat$W_b)

dat$beta<-c(66,558,280540,100, 113, 381,22,100)

dat$sigma<-c(1.3,2.1,3.2,0.5,1.6,1.9,1.5,0.5) ## nephrops standard from LemaRns ask shellfish 
dat$k_vb<-dat$k
dat$erepro<-0.2
dat$r_max<-c(8e9, 2e12, 1e12, 1e12, 4e14,4e10,5e11, 1e12)

dat1<-dat[,c("species","w_inf","w_mat", "beta", "sigma", "k_vb", "erepro", "r_max")]

irish_tau <- read.csv("Data/IS7_tau.csv")
colnames(irish_tau)[4] <- "Nephrops_W"
irish_tau <- cbind(rbind(irish_tau,irish_tau[4,]),0)
colnames(irish_tau)[8] <- "Nephrops_E"
sp_ol <- as.matrix(read.csv("Data/irish_SpOverlap.csv")[,-1])
irish_tau <- t(irish_tau * sp_ol)

parms<-MizerParams(dat1,irish_tau , kappa=1e11)
g_sel <- read.csv("Data/Fish_IS7_logistic.txt",skip = 1,header = F)
g_sel

is_L50 <- g_sel[,3];is_L50[8] <- is_L50[4]
is_eta <- g_sel[,2];is_eta[8] <- is_eta[4]

is_L25<- (-log(3)/is_eta)+ is_L50
is_L25


dat1$a<- dat$W_a
dat1$b<- dat$W_b
dat1$l50<-is_L50
dat1$l25<- is_L25
dat1$sel_func<-"sigmoid_length"
dat1

parms<-MizerParams(dat1,irish_tau , kappa=1e11)

sim <- project(parms)


plot(sim)
## building fishing dataframe
delta_t <- 0.1 #timestep  of the numerical integration 
effort<-array(0, dim=c((100+47)*1/(delta_t) + 1,8)) #effort matrix 
dimnames(effort)[2] <- list(parms@species_params$species) #labeling effort matrix by gear 
dimnames(effort)[1] <- list(seq(0,147,delta_t))# label time step 
sim<- project(parms,effort=effort,dt=delta_t)

### adding in fishing values from the assessment

fs <- read.csv("Data/assessed_fs.csv") ## Fs from the assessment
new_fs <- apply(fs[,-1],2,function(x){rep(x,each=1/(delta_t))}) ## sort them to put in the effort matrix
effort[-(1:(100 * 1/ delta_t + 1)),] <- new_fs ## put them in the effort matrix (not including the burn in)

sim<- project(parms,effort=effort,dt=delta_t) ## run the model
plot(sim) ## plot the results

## The model is ready to fit:
## 1) Change the effort for the first 100 years of the simulation (this is spin-up F and there is one per species)
spinUpF<- c(0.1,0.2,0.1,0.14,0.13,0.12,0.5,0.9)

effort[(1:(100 * 1/ delta_t + 1)),] <-matrix(spinUpF, 100*1/delta_t + 1, 8, byrow = TRUE)

## 2) Import the data from "Data/catches.csv"
catch<-read.csv("Data/catches.csv")
getYield(sim)
sim<- project(parms,effort=effort,dt=delta_t)

sim_catch<-log(getYield(sim)[-(1:(111)),]) - log(1e6) ### log(1e6) is converting grams to tonnes 

plot.ts(sim_catch)
plot.ts(log(catch))
## 3) Create a function that takes kappa, r_max and spin-up F and returns the modelled yield for the last 37 years of the simulation. 


inv_logit <- function(x){
  return(1 / (1 + exp(-x)))
}

cancatchthis_eps <-function(theta, effort, delta_t=0.1, dat1, irish_tau){ 
  dat1$r_max<- exp(theta[1:8])
  dat1$erepro <- exp(theta[9:16])
  #dat1$erepro <- inv_logit(theta[9:16])
  kappa <- exp(theta[17])
  spinUpF <- theta[18:25]
  parms<-MizerParams(dat1,irish_tau , kappa=kappa)
  effort[(1:(100 * 1/ delta_t+1)),] <-matrix(spinUpF, 100*1/delta_t + 1, 8, byrow = TRUE)
  sim<- project(parms,effort=effort,dt=delta_t)
  sim_catch<-log(getYield(sim)[-(1:(111)),]) - log(1e6)
  return(sim_catch)
}

system.time(cancatchthis_eps(theta = c(log(dat1$r_max),log(dat1$erepro), log(1e11),spinUpF), effort=effort, dat1 = dat1, irish_tau = irish_tau)) 
