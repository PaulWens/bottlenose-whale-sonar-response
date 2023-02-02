# Response Intensity modelling based on DTAG data

source("RI_functions.R") # fitting functions
library(bbmle)
library(stats4)

################ DATA

  load("DTAG-data.Rdata")
  
  # Exposures for ha13_176a, ha15_171a, ha15_179b, ha16_170a
  exp.RL <- c(151, 99, 128, 128) - 79 # minus hearing threshold from Pacini et al measured at 5 kHz

  #distances (km)
  exp.d <-c(4.4,0.02,0.8,16.8)
  
  #time of each one (index in table)
  exp.t <- which(tab$PB==1)
  
  # exposure identifier for each time step (0 for no exposure)
  exp.id <- rep(0, length(tab$id))
  for(j in 1:length(exp.t)) {
    exp.id[tab$id==tab$id[exp.t[j]] & tab$htot>=tab$htot[exp.t[j]]] <- j
  }
    
################ Model fitting

  # Define response variable (LL function takes these from workspace rather than as arguments)
  Mdist<-list(all=tab$dist)
  mdist <- tab$dist
  mdist_orig <- mdist
  T<-length(mdist) # total number of dives

  # beta1: b0 - RI intercept (expected MD when there is no exposure)
  # beta2: b1 - scales the initial response intensity as a function of received level
  # beta3: b2 - decay rate of the response
  # beta4: b3 - scales the RI as a function of range
  # k_gamma: shape parameter of gamma distribution

  ipv <- list(k_gamma=0.1,beta1=0.7,beta2=0.01,beta3=0.1,beta4=0.1)#initial values in natural space
  wp <- pn2pw(ipv) # initial values in working space
  
### AIC model selection
  MLM1<-list(all.full = mle2(LL, start = list(k_gamma=wp[1],beta1=wp[2],beta2=wp[3],beta3=wp[4],beta4=wp[5]),method = "L-BFGS-B")) 
  exp.RL <- c(1, 1, 1, 1) 
  MLM2<-list(all.full = mle2(LL, start = list(k_gamma=wp[1],beta1=wp[2],beta2=wp[3],beta3=wp[4],beta4=wp[5]),method = "L-BFGS-B")) 
  exp.RL <- c(151, 99, 128, 128) - 79
  MLM3<-list(all.full = mle2(LL, start = list(k_gamma=wp[1],beta1=wp[2],beta2=wp[3],beta4=wp[5]),method = "L-BFGS-B")) 
  MLM4<-list(all.full = mle2(LL, start = list(k_gamma=wp[1],beta1=wp[2],beta2=wp[3],beta3=wp[4]),method = "L-BFGS-B")) 
  exp.RL <- c(1, 1, 1, 1)
  MLM5<-list(all.full = mle2(LL, start = list(k_gamma=wp[1],beta1=wp[2],beta2=wp[3],beta4=wp[5]),method = "L-BFGS-B")) 
  MLM6<-list(all.full = mle2(LL, start = list(k_gamma=wp[1],beta1=wp[2],beta2=wp[3],beta3=wp[4]),method = "L-BFGS-B")) 
  exp.RL <- c(151, 99, 128, 128) - 79
  MLM7<-list(all.full = mle2(LL, start = list(k_gamma=wp[1],beta1=wp[2],beta2=wp[3]),method = "L-BFGS-B")) 
  MLM8<-list(all.full = mle2(LL, start = list(k_gamma=wp[1],beta1=wp[2]),method = "L-BFGS-B")) 
  
  getAIC(MLM1$all.full) # full model
  getAIC(MLM2$all.full) # intercept, time decay, range
  getAIC(MLM3$all.full) # intercept, RL, range
  getAIC(MLM4$all.full) # intercept, RL, decay
  getAIC(MLM5$all.full) # intercept, range
  getAIC(MLM6$all.full) # intercept, time decay
  getAIC(MLM7$all.full) # intercept, RL
  getAIC(MLM8$all.full) # intercept 

  # Estimates for best model
  bestModel <- MLM7
  print(bestModel$all.full)
  summary(bestModel$all.full)
  nparams <- pw2pn(as.numeric(coef(bestModel$all.full)),names(coef(bestModel$all.full)))
  nparams.l <- getHessianCI(bestModel$all.full)$lower
  nparams.u <- getHessianCI(bestModel$all.full)$upper

  # Plot observed vs fitted values
  D.maxL <- list(all=PredD(bestModel$all.full, mdist, exp.t=exp.t, exp.RL=exp.RL,exp.d=exp.d))
  plot(Mdist$all, 
       ylim=c(min(c(Mdist$all, D.maxL$all)), max(c(Mdist$all, D.maxL$all))), 
       type="b", pch=1+15*(exp.id>0))
  points(exp.t, Mdist$all[exp.t], col="red", pch=16)
  lines(D.maxL$all, col="darkgrey", lwd=2)
  
