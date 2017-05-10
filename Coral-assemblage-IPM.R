### IPM coral assemblage Dynamics ###
#DOI:https://doi.org/10.5281/zenodo.573772

#by Mohsen Kayal, May 2015
#mohsen.kayal@gmail.com

#this code uses empirical data on individual coral dynamics
#to estimate coral demographic performance in recruitment, growth & survival
#which is subsequently used to predict coral assemblage dynamics in a
#multi-species, open-population Integral Projection Model
#where coral recruitment is density-dependent
#and survival and growth are size-dependent

#this code was developed in application to coral dynamics data from Moorea
#where coral communities are dominated by an assemblage of 3 genera
#Acropora, Pocillopora & Porites

#for study design and approach, refer to Kayal et al. paper entitled
#"Predicting coral community dynamics and reef ecosystem resilience using population models"
#see Kayal et al. 2015 [Ecological Complexity] for quantitative approach to coral demography
#see Eller et al. 2016 [Springer] for IPM methodology & terminology




### 1. Estimate demographic parameters ###


# PREPARE DATA

Data <- read.table("DATA-FILE-PATH",header=T) #open data table

#Define sp. groups
DataAcr <- Data[which(SData$Taxa=="Acr"),] #Acr-opora
DataPoc <- Data[which(SData$Taxa=="Poc"),] #Poc-illopora
DataPor <- Data[which(SData$Taxa=="Por"),] #Por-ites
DataAss <- rbind(DataAcr, DataPoc, DataPor) #Ass-emblage

#Define variables for each group
#Acr
zAcr <- DataAcr$logSi #log initial size (cm2)
SurvAcr <- DataAcr$Survival #survival (binary)
z1Acr <- DataAcr$logSf #log final size (cm2)
#Poc
zPoc <- DataPoc$logSi
SurvPoc <- DataPoc$Survival
z1Poc <- DataPoc$logSf
#Por
zPor <- DataPor$logSi
SurvPor <- DataPor$Survival
z1Por <- DataPor$logSf



# CALCULATE SIZE-DEPENDENT SURVIVAL & GROWTH, AND DENSITY-DEPENDENT RECRUITMENT

#use a GLM function for survival data
mod.SurvAcr <- "GLM-FUNC"(y=SurvAcr, x=zAcr, random="RANDOM-EFFECTS", family=binomial, data=DataAcr)
mod.SurvPoc <- "GLM-FUNC"(y=SurvPoc, x=zPoc, random="RANDOM-EFFECTS", family=binomial, data=DataPoc)
mod.SurvPor <- "GLM-FUNC"(y=SurvPor, x=zPor, random="RANDOM-EFFECTS", family=binomial, data=DataPor)

#use a GLM function for growth data
mod.GrowAcr <- "GLM-FUNC"(y=z1Acr, x=zAcr, random="RANDOM-EFFECTS", family=gaussian, data=DataAcr)
mod.GrowPoc <- "GLM-FUNC"(y=z1Poc, x=zPoc, random="RANDOM-EFFECTS", family=gaussian, data=DataPoc)
mod.GrowPor <- "GLM-FUNC"(y=z1Por, x=zPor, random="RANDOM-EFFECTS", family=gaussian, data=DataPor)

#use a GLM function for recruitment data
#(not used in our paper as density-dependent recruitment was quantified in a previous study)
mod.RecrAcr <- "GLM-FUNC"(y=recr.densAcr, x=pop.densAcr, ...)
mod.RecrPoc <- "GLM-FUNC"(y=recr.densPoc, x=pop.densPoc, ...)
mod.RecrPor <- "GLM-FUNC"(y=recr.densPor, x=pop.densPor, ...)




# STORE DEMOGRAPHIC PARAMETERS

m.parAcr <- c(
  #size-dependent survival
  surv = coef(mod.SurvAcr), #GLM intercept & slope
  #size-dependent growth
  grow = coef(mod.GrowAcr), #GLM intercept & slope
  grow.sd = summary(mod.GrowAcr)$sigma, #GLM residual standard deviation
  #density-dependent recruitment
  recr = coef(mod.RecrAcr), #GLM intercept & slope
  #recruit size
  recr.zmean = "MEAN-SIZE", #can be predetermined or estimated from empirical data 
  recr.zsd = "SD-SIZE") #can be predetermined or estimated from empirical data

m.parPoc <- c(
  #size-dependent survival
  surv = coef(mod.SurvPoc),
  #size-dependent growth
  grow = coef(mod.GrowPoc),
  grow.sd = summary(mod.GrowPoc)$sigma,
  #density-dependent recruitment
  recr = coef(mod.RecrPoc),
  #recruit size
  recr.zmean = "MEAN-SIZE",
  recr.zsd = "SD-SIZE")

m.parPor <- c(
  #size-dependent survival
  surv = coef(mod.SurvPor),
  #size-dependent growth
  grow = coef(mod.GrowPor),
  grow.sd = summary(mod.GrowPor)$sigma,
  #density-dependent recruitment
  recr = coef(mod.RecrPor),
  #recruit size
  recr.zmean = "MEAN-SIZE",
  recr.zsd = "SD-SIZE")






### 2. Build demographic functions ###


# BUILD IPM FUNCTIONS

#function to calculate survival probability from size
S_zAcr <- function(zAcr, m.parAcr) {
  linear.pAcr <- m.parAcr["surv.Intercept"] + m.parAcr["surv.Slope"] * zAcr #linear predictor
  pAcr <- 1/(1+exp(-linear.pAcr)) #back-transformed probability
  return(pAcr)}
S_zPoc <- function(zPoc, m.parPoc) {
  linear.pPoc <- m.parPoc["surv.Intercept"] + m.parPoc["surv.Slope"] * zPoc
  pPoc <- 1/(1+exp(-linear.pPoc))
  return(pPoc)}
S_zPor <- function(zPor, m.parPor) {
  linear.pPor <- m.parPor["surv.Intercept"] + m.parPor["surv.Slope"] * zPor
  pPor <- 1/(1+exp(-linear.pPor))
  return(pPor)}

#function to calculate probability density function of future size from size
G_z1zAcr <- function(z1Acr, zAcr, m.parAcr) {
  muAcr <- m.parAcr["grow.Intercept"] + m.parAcr["grow.Slope"] * zAcr
  sigAcr <- m.parAcr["grow.sd"]
  p.den.growAcr <- dnorm(z1Acr, mean = muAcr, sd = sigAcr)
  return(p.den.growAcr)}
G_z1zPoc <- function(z1Poc, zPoc, m.parPoc) {
  muPoc <- m.parPoc["grow.Intercept"] + m.parPoc["grow.Slope"] * zPoc
  sigPoc <- m.parPoc["grow.sd"]
  p.den.growPoc <- dnorm(z1Poc, mean = muPoc, sd = sigPoc)
  return(p.den.growPoc)}
G_z1zPor <- function(z1Por, zPor, m.parPor) {
  muPor <- m.parPor["grow.Intercept"] + m.parPor["grow.Slope"] * zPor
  sigPor <- m.parPor["grow.sd"]
  p.den.growPor <- dnorm(z1Por, mean = muPor, sd = sigPor)
  return(p.den.growPor)}

#function to calculate number of new recruits from population density (on previous yr)
#used to estimate recruitment for Acr & Por in our study
R_zAcr <- function(z.distAcr, m.parAcr){
  linear.RAcr <- #linear predictor
    m.parAcr["recr.Intercept"]+m.parAcr["recr.Slope"]*prev.dAcr(z.distAcr, meshpts)
  RAcr <- exp(linear.RAcr) #back-transformed density
  return(RAcr) }
R_zPor <- function(z.distPor, m.parPor){
  linear.RPor <- 
    m.parPor["recr.Intercept"]+m.parPor["recr.Slope"]*prev.dPor(z.distPor, meshpts)
  RPor <- exp(linear.RPor)
  return(RPor) }

#function to calculate number of new recruits from assemblage surface (on previous yr)
#used to estimate recruitment for Poc in our study
R_zPoc <- function(z.distAss, m.parPoc){
  linear.RPoc <-
    m.parPoc["recr.Intercep"]+m.parPoc["recr.Slope"]*prev.surfAss(z.distAss, meshpts)
  RPoc <- exp(linear.RPoc)
  return(RPoc) }

#function to calculate population & assemblage density (on previous yr)
i=1 #define iterative variable
prev.dAcr <- prev.dPoc <- prev.dPor <- prev.dAss <- list() #create vector
prev.densityAcr <- function (z.distAcr, meshpts) {
  if(i==1) { prev.dAcr[[i]]=0 } else { #density=0 prior 1st yr
    prev.dAcr[[i]] <- sum(z.distAcr[[i-1]]*h.mesh) }
  return(prev.dAcr[[i]]) }
prev.densityPoc <- function (z.distPoc, meshpts) {
  if(i==1) { prev.dPoc[[i]]=0 } else {
    prev.dPoc[[i]] <- sum(z.distPoc[[i-1]]*h.mesh) }
  return(prev.dPoc[[i]]) }
prev.densityPor <- function (z.distPor, meshpts) {
  if(i==1) { prev.dPor[[i]]=0 } else {
    prev.dPor[[i]] <- sum(z.distPor[[i-1]]*h.mesh) }
  return(prev.dPor[[i]]) }
prev.densityAss <- function (z.distAss, meshpts) {
  if(i==1) { prev.dAss[[i]]=0 } else {
    prev.dAss[[i]] <- sum(z.distAss[[i-1]]*h.mesh) }
  return(prev.dAss[[i]]) }

#function to calculate population & assemblage surface (on previous yr)
i=1 #define iterative variable
prev.surfAcr <- prev.surfPoc <- prev.surfPor <- prev.surfAss <- list() #create vector
prev.surfaceAcr <- function (z.distAcr, meshpts) {
  if(i==1) { prev.surfAcr[[i]]=0 } else { #surface=0 prior 1st yr
    prev.surfAcr[[i]] <- sum(z.distAcr[[i-1]]*h.mesh*meshpts) }
  return(prev.surfAcr[[i]]) }
prev.surfacePoc <- function (z.distPoc, meshpts) {
  if(i==1) { prev.surfPoc[[i]]=0 } else {
    prev.surfPoc[[i]] <- sum(z.distPoc[[i-1]]*h.mesh*meshpts) }
  return(prev.surfPoc[[i]]) }
prev.surfacePor <- function (z.distPor, meshpts) {
  if(i==1) { prev.surfPor[[i]]=0 } else {
    prev.surfPor[[i]] <- sum(z.distPor[[i-1]]*h.mesh*meshpts) }
  return(prev.surfPor[[i]]) }
prev.surfaceAss <- function (z.distAss, meshpts) {
  if(i==1) { prev.surfAss[[i]]=0 } else {
    prev.surfAss[[i]] <- sum(z.distAss[[i-1]]*h.mesh*meshpts) }
  return(prev.surfAss[[i]]) }

#calculate recruit size
c_z1zAcr <- function(meshpts, m.parAcr) {
  mu <- m.parAcr["recr.zmean"]          #mean
  sig <- m.parAcr["recr.zsd"]           #sd
  p.den.recrzAcr <- dnorm(meshpts, mean = mu, sd = sig) #probability density function
  return(p.den.recrzAcr) }
c_z1zPoc <- function(meshpts, m.parPoc) {
  mu <- m.parPoc["recr.zmean"]
  sig <- m.parPoc["recr.zsd"]
  p.den.recrzPoc <- dnorm(meshpts, mean = mu, sd = sig)
  return(p.den.recrzPoc) }
c_z1zPor <- function(meshpts, m.parPor) {
  mu <- m.parPor["recr.zmean"]
  sig <- m.parPor["recr.zsd"]
  p.den.recrzPor <- dnorm(meshpts, mean = mu, sd = sig)
  return(p.den.recrzPor) }




# BUILD IPM KERNELS

#Define survival-growth kernel
P_z1zAcr <- function (z1Acr, zAcr, m.parAcr) {
  return(S_zAcr(zAcr, m.parAcr) * G_z1zAcr(z1Acr, zAcr, m.parAcr)) }
P_z1zPoc <- function (z1Poc, zPoc, m.parPoc) {
  return(S_zPoc(zPoc, m.parPoc) * G_z1zPoc(z1Poc, zPoc, m.parPoc)) }
P_z1zPor <- function (z1Por, zPor, m.parPor) {
  return(S_zPor(zPor, m.parPor) * G_z1zPor(z1Por, zPor, m.parPor)) }

#Define reproduction kernel
F_z1zAcr <- function (z.distAcr, meshpts, m.parAcr) { #(varies with pop. density)
  return( R_zAcr(z.distAcr, m.parAcr) * c_z1zAcr(meshpts, m.parAcr) ) }
F_z1zPoc <- function (z.distAss, meshpts, m.parPoc) { #(varies with ass. surface)
  return( R_zPoc(z.distAss, m.parPoc) * c_z1zPoc(meshpts, m.parPoc) ) }
F_z1zPor <- function (z.distPor, meshpts, m.parPor) { #(varies with pop. density)
  return( R_zPor(z.distPor, m.parPor) * c_z1zPor(meshpts, m.parPor) ) }

#Build discretized kernel using mesh-points
mk_KAcr <- function(m, m.parAcr, L, U) {
  #create mesh points
  L <- 0; U <- 7 #Lower & Upper size limits)
  m <- 100 #number of mesh-points
  h.mesh <- (U-L)/m #distance between mesh-points
  meshpts <- L + (1:m)*h.mesh - h.mesh/2
  #compute iteration matrix
  PAcr <- h.mesh * (outer(meshpts, meshpts, P_z1zAcr, m.parAcr = m.parAcr))
  FAcr <- h.mesh * (outer(meshpts, meshpts, F_z1zAcr, m.parAcr = m.parAcr))
  KAcr <- PAcr + FAcr #full IPM kernell
  return(list( KAcr=KAcr, meshpts=meshpts, PAcr=PAcr, FAcr=FAcr, h.mesh=h.mesh )) }
mk_KPoc <- function(m, m.parPoc, L, U) {
  #create mesh points
  L <- 0; U <- 7
  m <- 100
  h.mesh <- (U-L)/m
  meshpts <- L + (1:m)*h.mesh - h.mesh/2
  #compute iteration matrix
  PPoc <- h.mesh * (outer(meshpts, meshpts, P_z1zPoc, m.parPoc = m.parPoc))
  FPoc <- h.mesh * (outer(meshpts, meshpts, F_z1zPoc, m.parPoc = m.parPoc))
  KPoc <- PPoc + FPoc
  return(list( KPoc=KPoc, meshpts=meshpts, PPoc=PPoc, FPoc=FPoc, h.mesh=h.mesh )) }
mk_KPor <- function(m, m.parPor, L, U) {
  #create mesh points
  L <- 0; U <- 7
  m <- 100
  h.mesh <- (U-L)/m
  meshpts <- L + (1:m)*h.mesh - h.mesh/2
  #compute iteration matrix
  PPor <- h.mesh * (outer(meshpts, meshpts, P_z1zPor, m.parPor = m.parPor))
  FPor <- h.mesh * (outer(meshpts, meshpts, F_z1zPor, m.parPor = m.parPor))
  KPor <- PPor + FPor
  return(list( KPor=KPor, meshpts=meshpts, PPor=PPor, FPor=FPor, h.mesh=h.mesh )) }



# COMBINE IPM KERNELS

nBigMatrix <- 100

## set lower and upper limits for the size range
min.size <- 0.01
max.size <- 7

#IPM kernels
IPM.Acr <- mk_KAcr(nBigMatrix,m.parAcr,min.size,max.size)
IPM.Poc <- mk_KPoc(nBigMatrix,m.parPoc,min.size,max.size)
IPM.Por <- mk_KPor(nBigMatrix,m.parPor,min.size,max.size)






### 3. Run IPM simulations ###


# PREDICT REEF COLONISATION BY A CORAL ASSEMBLAGE DRIVEN BY LARVAL RECRUITMENT

#extract mesh point information from one of the kernels
meshpts <- IPM.Acr$meshpts
h.mesh <- IPM.Acr$h.mesh

#create size probability-density-distribution vectors
z.distAcr <- z.distPoc <- z.distPor <- z.distAss <- list()

i=1 #reset iterative variable

#on 1st yr, assemblage is made of new recruits from the 3 populations
  z.distAcr[[1]] <- exp(m.parAcr["recr.Intercept"]) * c_z1zAcr(meshpts,m.parAcr)
  z.distPoc[[1]] <- exp(m.parPoc["recr.Intercept"]) * c_z1zPoc(meshpts,m.parPoc)
  z.distPor[[1]] <- exp(m.parPor["recr.Intercept"]) * c_z1zPor(meshpts,m.parPor)
  z.distAss[[1]] <- z.distAcr[[1]]+z.distPoc[[1]]+z.distPor[[1]]

#for following 5 yrs, size-dependent survival & growth + density-dependent recruitment
  for (i in 2:6) {
    #Acr
    z.distAcr[[i]] <- IPM.Acr$PAcr %*% z.distAcr[[i-1]] + 
      F_z1zAcr(z.distAcr,meshpts,m.parAcr)
    #Poc
    z.distPoc[[i]] <- IPM.Poc$PPoc %*% z.distPoc[[i-1]] + 
      F_z1zPoc(z.distAss, meshpts, m.parPoc)
    #Por
    z.distPor[[i]] <- IPM.Por$PPor %*% z.distPor[[i-1]] + 
      F_z1zPor(z.distPor,meshpts,m.parPor)
    #Ass
    z.distAss[[i]] <- z.distAcr[[i]]+z.distPoc[[i]]+z.distPor[[i]] }

#plot changes in size probability-density-distributions across yrs
#Acr
plot(z.distAcr[[1]]~meshpts,col="blue",type="o",pch="1",cex=0.5,
     ylim=c(0,max(unlist(z.distAcr))),xlab="Coral log-size",
     ylab="Probability density distribution by year")
lines(z.distAcr[[2]]~meshpts, col="blue", type="o", pch="2", cex=0.5)
lines(z.distAcr[[3]]~meshpts, col="blue", type="o", pch="3", cex=0.5)
lines(z.distAcr[[4]]~meshpts, col="blue", type="o", pch="4", cex=0.5)
lines(z.distAcr[[5]]~meshpts, col="blue", type="o", pch="5", cex=0.5)
lines(z.distAcr[[6]]~meshpts, col="blue", type="o", pch="6", cex=0.5)
#Poc
plot(z.distPoc[[1]]~meshpts,col="red",type="o",pch="1",cex=0.5,
     ylim=c(0,max(unlist(z.dist.by.yrPoc))),xlab="Coral log-size",
     ylab="Probability density distribution by year")
lines(z.distPoc[[2]]~meshpts, col="red", type="o", pch="2", cex=0.5)
lines(z.distPoc[[3]]~meshpts, col="red", type="o", pch="3", cex=0.5)
lines(z.distPoc[[4]]~meshpts, col="red", type="o", pch="4", cex=0.5)
lines(z.distPoc[[5]]~meshpts, col="red", type="o", pch="5", cex=0.5)
lines(z.distPoc[[6]]~meshpts, col="red", type="o", pch="6", cex=0.5)
#Por
plot(z.distPor[[1]]~meshpts,col="green",type="o",pch="1",cex=0.5,
     ylim=c(0,max(unlist(z.dist.by.yrPor))),xlab="Coral log-size",
     ylab="Probability density distribution by year")
lines(z.distPor[[2]]~meshpts, col="green", type="o", pch="2", cex=0.5)
lines(z.distPor[[3]]~meshpts, col="green", type="o", pch="3", cex=0.5)
lines(z.distPor[[4]]~meshpts, col="green", type="o", pch="4", cex=0.5)
lines(z.distPor[[5]]~meshpts, col="green", type="o", pch="5", cex=0.5)
lines(z.distPor[[6]]~meshpts, col="green", type="o", pch="6", cex=0.5)
#Ass
plot(z.distAss[[1]]~meshpts,col="black",type="o",pch="1",cex=0.5,
     ylim=c(0,max(unlist(z.dist.by.yrAss))),xlab="Coral log-size",
     ylab="Probability density distribution by year")
lines(z.distAss[[2]]~meshpts, col="black", type="o", pch="2", cex=0.5)
lines(z.distAss[[3]]~meshpts, col="black", type="o", pch="3", cex=0.5)
lines(z.distAss[[4]]~meshpts, col="black", type="o", pch="4", cex=0.5)
lines(z.distAss[[5]]~meshpts, col="black", type="o", pch="5", cex=0.5)
lines(z.distAss[[6]]~meshpts, col="black", type="o", pch="6", cex=0.5)







### 4. Estimate elasticity of population trajectories ###


#  CALCULATE PROPORTIONAL CHANGE IN POPULATION TRAJECTORY

#create function that calculates proportional change in pop. surface after 5 yrs
#for a population whose demographic parameters are stored in a vector x
#(here with example of Acr dynamics) 
surf.changeAcr <- function(x) {
  IPM.Acr <- mk_KAcr(nBigMatrix,x,min.size,max.size) #build kernel
  
  #generate yr 1
  z.distAcr <- list() 
  z.distAcr[[1]] <- exp(x["recr.Intercept"]) * c_z1zAcr(meshpts,x)
  
  for (i in 2:6) { #generate consecutive 5 yrs
    
    #density (on previous yr) function
    if(i==1) { prev.dAcr[[i]]=0 } else {
      prev.dAcr[[i]] <- sum(z.distAcr[[i-1]]*h.mesh) }
    
    #density-dependent recruitment function
    linear.RAcr <- list()
    linear.RAcr[[i]] <- x["recr.Intercept"]+x["recr.Slope"]*prev.dAcr[[i]]
    RAcr[[i]] <- exp(linear.RAcr[[i]])

    #run simulation
    z.distAcr[[i]] <- 
      IPM.Acr$PAcr %*% z.distAcr[[i-1]] + RAcr[[i]]* c_z1zAcr(meshpts,x) }
  
  #calculate variation in pop. surface
  proj.surfAcr <- list()
  for(i in 1:6) proj.surfAcr[[i]]<-sum(z.distAcr[[i]]*h.mesh*meshpts) #pop. surface
  surf.ratioAcr<-proj.surfAcr[[6]]/proj.surfAcr[[1]] #ratio in pop. surface
  return(surf.ratioAcr) }

#apply function to the demographic parameters estimated for Acr
surf.changeAcr(m.parAcr)



# CALCULATE ELASTICITY OF POPULATION TRAJECTORY TO MODEL PARAMETERS

#Calculate sensitivity of population surface to demographic parameters
library(numDeriv)
sensitivity.surf.changeAcr <- grad(surf.changeAcr, m.parAcr)

#Calculate elasticity of population surface to demographic parameters
elasticity.surf.changeAcr <- abs(sensitivity.surf.changeAcr*m.parAcr/surf.changeAcr(m.parAcr))


###End of IPM coral assemblage dynamics code
