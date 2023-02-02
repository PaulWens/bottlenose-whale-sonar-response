# Functions by S DeRuiter and S Isojunno for fitting Response Intensity models to MD data

#######################################
#  Calculate expected values of model
#######################################
#Return the expected values of Distance as a function of time for a model with a specified parameter set
#t is a scalar time (time block index)
#exp.RL is a vector of RLs during exposed dives (only)
#exp.t is a vector of time blocks during which exposures occurred
#beta1,2,3,4 are model parameters (some may be omitted)

#note that model parameter names are slightly different in R-code:
#   TEXT NAME     R-code NAME
#   omega         k_gamma
#   beta_0        beta1
#   beta_1        beta2
#   beta_2        beta3
#   beta_3        beta4

E.D <-function(t,exp.RL,exp.t,exp.d, beta1,beta2,beta3,beta4){#input params can be any subset of beta1, beta2, beta3, beta4. if a parameter is not input, or is NA, the function assumes that the corresponding term should be omitted from the model.
  #input check - if any inputs are missing, set to NA
  beta1 <- ifelse(missing(beta1),NA,beta1)
  beta2 <- ifelse(missing(beta2),NA,beta2)
  beta3 <- ifelse(missing(beta3),NA,beta3)
  beta4 <- ifelse(missing(beta4),NA,beta4)
  #calculations  
  res   <- ifelse(is.na(beta1),0,beta1) #baseline average distance is beta1 if beta1 exists, 0 otherwise
  if(exp.id[t]>0) {
    e0    <- matrix(c(exp.RL*beta2, exp(beta3*(exp.t-t)), (1/(1+exp.d*beta4))), 
                    nrow=3, ncol=length(exp.t), byrow=TRUE)
    effect<- apply(e0,MARGIN=2, FUN=prod, na.rm=TRUE)
    effect[effect==1] <- 0  ## ADDED BASED ON COMMENT OF SD 22-sep-2016 (prod() returns 1 if vector is all NA and na.rm=TRUE)
    res<-res+sum(effect[exp.id[t]])
  }
  return(res)
}

############################################################
# Calculate minus log-likelihood of a model given the data
############################################################

LL<-function(k_gamma,beta1,beta2=-Inf,beta3=0,beta4=-Inf){
  nparams <- pw2pn(c(k_gamma,beta1,beta2,beta3,beta4))
  D<-vapply(X=seq(from=1,to=length(mdist),by=1) , 
            FUN=E.D, FUN.VALUE=1, exp.RL=exp.RL,exp.t=exp.t, exp.d=exp.d, 
            beta1=nparams$beta1, beta2=nparams$beta2, beta3=nparams$beta3, beta4=nparams$beta4)
  k_gamma[k_gamma==0] <- 1e-9 #to prevent log(0), /0
  D[D==0] <- 1e-9#to prevent log(0), /0
  s<-nparams$k_gamma/D
  a<-(D^2)/nparams$k_gamma
  a[a==0] <- a[a==0]+1e-9 #to prevent "value out of range in lgamma"
  a[ a<0 & is.integer(a)] <- a[ a<0 & is.integer(a)]+1e-9 #to prevent "value out of range in lgamma"
  LLike <- sum((mdist/s)-(a-1)*log(mdist)+lgamma(a)+a*log(s))
  return(LLike)
}

###############################################################################
# Get fitted values for a given model (without calculating the likelihood)
###############################################################################
PredD <- function(FittedModel, mdist, exp.t, exp.RL, exp.d)#get predicted distance values (fitted model)
{
  nparams <- pw2pn(as.numeric(coef(FittedModel)),names(coef(FittedModel)))
  if(length(nparams$beta3)==0) {nparams$beta3 <- 0}
  if(length(nparams$beta4)==0) {nparams$beta4 <- 0}
  D <- vapply(X=seq(from=1,to=length(mdist),by=1) , FUN=E.D, FUN.VALUE=1, exp.RL=exp.RL,exp.t=exp.t, exp.d=exp.d, 
              beta1=nparams$beta1, beta2=nparams$beta2, beta3=nparams$beta3, beta4=nparams$beta4)
  return(D)
}

###############################################################################
# Predict values, given a set of parameters
###############################################################################
PredD2 <- function(predcoefs, mdist, exp.t, exp.RL, exp.d)#get predicted distance values (fitted model)
{
  if(!any(names(predcoefs)=="beta3")) {predcoefs$beta3 <- 0}
  D <- vapply(X=seq(from=1,to=length(mdist),by=1) , FUN=E.D, FUN.VALUE=1, 
              exp.RL=exp.RL,exp.t=exp.t, exp.d=exp.d, 
              beta1=predcoefs$beta1, beta2=predcoefs$beta2, beta3=predcoefs$beta3, beta4=predcoefs$beta4)
  return(D)
}

####################################################
# helper functions - transform parameters etc
####################################################
## function that transforms each of the (possibly constrained) parameters to 
#the real line (and converts list to a vector)
pn2pw <- function(nparams)
{
  lpar <- lapply(nparams,log)
  wparvect <- as.numeric(lpar)
  wparvect[names(lpar)=="beta3"] <- nparams["beta3"][[1]] # beta3 not transformed
  return(wparvect)
}

#inverse transform back to natural parameter space
# (and converts vector to a list)
pw2pn <- function(wparams,pnames=c("k_gamma", "beta1", "beta2", "beta3", "beta4")) 
{
  epar <- exp(wparams)
  epar[pnames=="beta3"] <- wparams[pnames=="beta3"] # beta3 not transformed
  lpar <- as.list(epar)
  names(lpar) <- pnames
  return(lpar)
}
#function to get parameters, in natural parameter space, from fitted model object
getnp <- function(FittedModel)
{
  ncoef <- pw2pn(as.numeric(coef(FittedModel)), names(coef(FittedModel))) #parameters in natural space
  return(ncoef)
}

#function to get parameter Std errors out in spite of the transformation of params
getSEs <- function(FittedModel,N)#N is sample size, here length(MDist)
{
  ncoef <- getnp(FittedModel)#parameters in natural space
  ses <- sqrt((diag(vcov(FittedModel)))/N)#these SEs are correct for parameters NOT transformed 
  #use delta method (first term of Taylor series approx) to get SEs for transformed parameters.
  ses[names(ses)=="k_gamma"] <- ses[names(ses)=="k_gamma"]*ncoef$k_gamma
  ses[names(ses)=="beta1"] <- ses[names(ses)=="beta1"]*ncoef$beta1
  ses[names(ses)=="beta2"] <- ses[names(ses)=="beta2"]*ncoef$beta2
  ses[names(ses)=="beta4"] <- ses[names(ses)=="beta4"]*ncoef$beta4
  return(ses)  
}

#function to get CIs in natural parameter space
getCIs <- function(FittedModel, level, parm=c("k_gamma", "beta1", "beta2", "beta3", "beta4"))#level is e.g. 0.95 for 95% CI
{
  ci0 <- confint(FittedModel, level=level, parm=parm)
  #transform k_gamma, beta1, beta2, beta4. if they are there.
  ci0[rownames(ci0)=="k_gamma"] <- exp(ci0[rownames(ci0)=="k_gamma"])
  ci0[rownames(ci0)=="beta1"] <- exp(ci0[rownames(ci0)=="beta1"])
  ci0[rownames(ci0)=="beta2"] <- exp(ci0[rownames(ci0)=="beta2"])
  ci0[rownames(ci0)=="beta4"] <- exp(ci0[rownames(ci0)=="beta4"])
  return(ci0)
}


getAIC <- function(FittedModel) {return(2*length(FittedModel@coef)+2*FittedModel@min)}

getHessianCI <- function(FittedModel, working=F) {
  wparams <- coef(FittedModel)
  pnames <- names(wparams)
  nparams <- pw2pn(as.numeric(wparams),pnames)
  ses.w <- sqrt(diag(solve(FittedModel@details$hessian)))
  ciw.l <- summary(FittedModel)@coef[,1]-qnorm(0.975)*ses.w
  ciw.u <- summary(FittedModel)@coef[,1]+qnorm(0.975)*ses.w
  cin.l <- pw2pn(as.numeric(ciw.l),pnames)
  cin.u <- pw2pn(as.numeric(ciw.u),pnames)
  if(working) {
    return(list(est=wparams, lower=ciw.l, upper=ciw.u))
  } else {
    return(list(est=nparams, lower=cin.l, upper=cin.u))
  }
}
