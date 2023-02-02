# Fit 3-state HMM models to ARGOS tracks filtered with crawl

Sys.setenv(tz="UTC")

# Load R Packages
library(crawl)
library(sp)
library(lubridate)
library(momentuHMM) # version 1.4.2

#### PREPARE DATA
load("sat-tag-data.RData") # load filtered track data
bwData <- prepData(data=crwOut) # preprocess data

# Find index of CEE start in each track
isexp <- rep(0,13)
for(id in 5:13){
  if(id<8){
    isexp[id] <- min(which(bwData$unique_posix >= tsexp15 & bwData$ID == id))
  }else{
    isexp[id] <- min(which(bwData$unique_posix >= tsexp16 & bwData$ID == id))
  }
}

# Add variable dist (source-whale distance at start of exposure)
poss15 <- data.frame(lat=71.01955, lon=-6.345034) # 2015-2
poss16 <- data.frame(lat=70.76058, lon=-6.096467) # 2016
poss <- rbind(poss15, poss16)
llcoord <- SpatialPoints(poss[,c(2,1)], proj4string=CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=29 ellps=WGS84"))
poss <- matrix(attr(utmcoord,"coords"), nrow=2)
row.names(poss) <- c("source15","source16")
tmp <- prepData(data=crwOut, centers=poss)
bwData$dist <- rep(0,length(bwData$ID))
for(id in 5:13){
  if(id<8){
    bwData$dist[bwData$ID==id] <- tmp$source15.dist[isexp[id]]
  }else{
    bwData$dist[bwData$ID==id] <- tmp$source16.dist[isexp[id]]
  }
}

# Add variable time_to_recovery
trecovmax = 8 # Maximum recovery time in hours (based on AIC)
timestep = 1 # Timestep in hours
bwData$time_to_recovery <- rep(0,length(bwData$ID))
trecov <- rev(seq(0,trecovmax,timestep)) # time to maximum recovery time
ind <- seq(0,length(trecov)-1)
for(id in c(5:13)){
  bwData$time_to_recovery[isexp[id]+ind] <- trecov
}

# Add variable SPLmax
RL = c(82,NA,NA,121,124,126,121,120,122) # SPLmax based on prop modelling (medians)
rs = tmp$source15.dist[isexp[c(6,7)]] 
RL[2:3] <- 185 - 20*log10(rs) - 3.6e-5*1^1.5*rs # SPLmax based on simple model
bwData$SPLmax <- rep(0,length(bwData$ID))
for(id in 5:13){
  bwData$SPLmax[bwData$ID==id] <- RL[id-4]
}  

# Remove step length and turn angle for segments without observed locations
tmp <- round(bwData$step,4)
ind <- diff(tmp)==0 # based on equal step lengths
bwData$step[ind] <- NaN
bwData$angle[ind] <- NaN

###### OPTIONAL: Remove individuals as a robustness check
#bwData <- bwData[!(bwData$ID==10|bwData$ID==12),]# #9=161588, 10=161590, 12=161592, 13=161593

###### Fit HMMs
stateNames = c("state 1: tortuous","state 2: low speed / directional","state 3: high speed / directional") # state labels

# distributions for observation processes
dist = list(step = "gamma", angle = "vm")

## Model 1 BASELINE MODEL. no covariates
nbStates <- 3 # number of states
step0 = rep(c(mean(bwData$step,na.rm=T), sd(bwData$step,na.rm=T)), nbStates)
#step0 = c(2500,2000,5000,3000,5000,3000) # some tuning needed when excluding IDs 9 and 13
angle0 = c(1,8,25) # low value = nondirectional 
m1 <- fitHMM(data = bwData, nbStates = nbStates, dist = dist, 
             Par0 = list(step=step0, angle=angle0),  
             retryFits=0, stateNames=stateNames,nlmPar=list(steptol=1e-6,print.level=0))

## Model 2 EFFECT OF EXPOSURE
fixbeta <- matrix(c(NA,NA,NA,NA,NA,0,NA,0,NA,0,NA,0), nrow=2) # only first row in TPM, from state 1 to 1,2,3
formula <- ~ time_to_recovery
Par0 <- getPar0(model=m1, nbStates=nbStates,
                estAngleMean=list(angle=F),
                circularAngleMean=list(angle=F), formula=formula)
m2 <- fitHMM(data = bwData, nbStates = nbStates, dist = dist, 
             Par0 = list(step=Par0$Par$step, angle=Par0$Par$angle),
             beta0 = Par0$beta, formula = formula,
             fixPar=list(beta=fixbeta),
             estAngleMean = list(angle=F), circularAngleMean = list(angle=F),
             retryFits=0, stateNames=stateNames, nlmPar=list(steptol=1e-6, print.level=0)) 

## Model 3 EFFECT OF RL
formula <- ~ time_to_recovery:SPLmax
Par0 <- getPar0(model=m1, nbStates=nbStates,
                estAngleMean=list(angle=F),
                circularAngleMean=list(angle=F), formula=formula)
m3 <- fitHMM(data = bwData, nbStates = nbStates, dist = dist, 
               Par0 = list(step=Par0$Par$step, angle=Par0$Par$angle),
               beta0 = Par0$beta, formula = formula,
               fixPar=list(beta=fixbeta),
               estAngleMean = list(angle=F), circularAngleMean = list(angle=F),
               retryFits=0, stateNames=stateNames,nlmPar=list(steptol=1e-6, print.level=0)) 

## Model 4 EFFECT OF DISTANCE
formula <- ~ time_to_recovery:dist
Par0 <- getPar0(model=m1, nbStates=nbStates,
                estAngleMean=list(angle=F),
                circularAngleMean=list(angle=F), formula=formula)
m4 <- fitHMM(data = bwData, nbStates = nbStates, dist = dist, 
               Par0 = list(step=Par0$Par$step, angle=Par0$Par$angle),
               beta0 = Par0$beta, formula = formula,
               fixPar=list(beta=fixbeta),
               estAngleMean = list(angle=F), circularAngleMean = list(angle=F),
               retryFits=0, stateNames=stateNames, nlmPar=list(steptol=1e-6, print.level=0))
AIC(m1, m2, m3, m4)

###### PLOT RESULTS
plot(m2,plotCI=T,plotTracks=T,ask=F,col=c("red","green","black"),cumul=F,covs=data.frame(time_to_recovery=max(bwData$time_to_recovery), SPLmax=129, dist=16100))

# plot pseudo-residuals for the steps
plotPR(m2,crash)

# Decode most likely state sequence
bwData$states <- viterbi(m2)

