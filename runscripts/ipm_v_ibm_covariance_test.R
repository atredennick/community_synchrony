##  Script to calculate IPM transition matrices from IBM population vector.
##  This is to compare the covariance matrix from the IBM and the IPM approx.

rm(list=ls()) # Clear the workspace

####
#### Load libraries
####
library(communitySynchrony)
library(synchrony)
library(boot)
library(mvtnorm)
library(msm)
library(MASS)
library(statmod)


####
####  Set global parameters
####
sppList=c("ARTR","HECO","POSE","PSSP")
bigM=c(50,75,50,75)       # Set matrix dimension for each species
maxSize=c(3000,202,260,225) # In cm^2: PSSP=225 HECO=202  POSE=260  ARTR=3000 
A=10000                # Area of 100cm x 100cm quadrat
Nspp=length(sppList)
tlimit=2
constant=T
Nyrs=22
NoOverlap.Inter=F


####
####  Read in regression parameters, subset for Idaho and format
####
do_site <- "Idaho"
Gpars_all <- readRDS("../results/growth_params_list.RDS")[[do_site]]
Spars_all <- readRDS("../results/surv_params_list.RDS")[[do_site]]
Rpars_all <- readRDS("../results/recruit_parameters.RDS")[[do_site]]

spp_list <- sppList <- names(Gpars_all)
Nyrs <- nrow(Gpars_all[[1]])
Nspp <- length(spp_list)

site_path <- paste("../data/", do_site, sep="")
Gpars <- format_growth_params(do_site = do_site, species_list = spp_list, 
                              Nyrs = Nyrs, Gdata_species = Gpars_all)
Spars <- format_survival_params(do_site = do_site, species_list = spp_list, 
                                Nyrs = Nyrs, Sdata_species = Spars_all)
Rpars <- format_recruitment_params(do_site = do_site, species_list = spp_list, 
                                   Nyrs = Nyrs, Rdata_species = Rpars_all,
                                   path_to_site_data = site_path)

if(constant==T){
  #turn off random year effects  
  Rpars$intcpt.yr=matrix(Rpars$intcpt.mu,Nyrs,Nspp,byrow=T)
  Gpars$intcpt.yr[]=0;Gpars$slope.yr[]=0
  Spars$intcpt.yr[]=0;Spars$slope.yr[]=0    
}


####
####  Vital rate functions
####
S=function(u,W,Spars,doYear,doSpp){
  mu=Spars$intcpt[doSpp]+Spars$intcpt.yr[doYear,doSpp]+
    (Spars$slope[doSpp]+Spars$slope.yr[doYear,doSpp])*u+
    W%*%(Spars$nb[doSpp,])
  return(inv.logit(mu))
}

G=function(v,u,W,Gpars,doYear,doSpp){
  mu=Gpars$intcpt[doSpp]+Gpars$intcpt.yr[doYear,doSpp]+(Gpars$slope[doSpp]+Gpars$slope.yr[doYear,doSpp])*u+
    W%*%(Gpars$nb[doSpp,])
  sigma2=Gpars$sigma2.a[doSpp]*exp(Gpars$sigma2.b[doSpp]*mu)
  out=dnorm(v,mu,sqrt(sigma2))
  out
}

#number of recruits per area produced 
# cover is stored in absolute area (cm^2)
get.rpa=function(Rpars,cover,doYear){
  # cover is in m^2 per m^2; convert to % scale:
  cover2=cover*100
  # calculate recruits
  Nspp=length(cover)
  mu=rep(NA,Nspp)
  for(i in 1:Nspp){
    mu[i]=cover2[i]*exp(Rpars$intcpt.yr[doYear,i]+sqrt(cover2)%*%Rpars$dd[i,]) 
  }
  if(sum(is.na(mu))>0) browser() # stop for errors
  rpa=mu/(cover*A)  # convert from number recruits to recruits per cm^2
  return(rpa)
}

# Fecundity function, expected number of recruits of size y produced by a size x individual
# The size distribution of recruits is on the log scale
f=function(v,u,Rpars,rpa,doSpp) { 
  nRecruits = rpa[doSpp]*exp(u)
  #probability of producing a seedling of size v
  tmp=dnorm(v,Rpars$sizeMean[doSpp],sqrt(Rpars$sizeVar[doSpp]))/(1-pnorm(-1.61,Rpars$sizeMean[doSpp],sqrt(Rpars$sizeVar[doSpp])))
  #number recruits of each size 
  f=nRecruits*tmp
  return(f)
}   



####
####  Read in IBM population vector list; keep iteration 50 as initial pop vector
####
nt.save <- readRDS("../results/nt_popvec_ibm.RDS")
nt <- list()
nt[[1]] <- nt.save[[2]][1:length(nt.save[[2]])]
nt[[1]][] <- 0
for(i in 1:length(nt.save)) nt[[i+1]] <- nt.save[[i]][1:length(nt.save[[i]])]


####
####  Simulation length, Matrix size and initial vectors
####
v=v.r=b.r=expv=Cr=WmatG=WmatS=list(length(sppList))
h=r.L=r.U=Ctot=numeric(length(sppList))
for(i in 1:Nspp){
  # minimum (0.9*minimum size from data) and maximum sizes (1.1*maximum size from data)
  L=log(0.2)
  U=log(maxSize[i])*1.1     
  # boundary points b and mesh points y. Note: b chops up the size interval (L-U) into bigM-equal-sized portions.
  b = L+c(0:bigM[i])*(U-L)/bigM[i] 
  # v calculates the middle of each n-equal-sized portion.
  v[[i]] = 0.5*(b[1:bigM[i]]+b[2:(bigM[i]+1)])
  # step size for midpoint rule. (see equations 4 and 5 in Ellner and Rees (2006) Am Nat.)
  h[i] = v[[i]][2]-v[[i]][1]  
  # variables for Wr approximation
  b.r[[i]]=sqrt(exp(b)/pi)
  v.r[[i]]=sqrt(exp(v[[i]])/pi)
  expv[[i]]=exp(v[[i]])
  r.L[i] = sqrt(exp(L)/pi)
  r.U[i] = sqrt(exp(U)/pi)
  WmatG[[i]]=matrix(NA,length(v.r[[i]]),Nspp)  # storage of size-specific W values for each focal species
  WmatS[[i]]=matrix(NA,length(v.r[[i]]),Nspp)
} # next species
tmp=range(v.r)
size.range=seq(tmp[1],tmp[2],length=50) # range across all possible sizes


####
####  Utility functions
####
get_pairs <- function(X, pop_vector){
  pairs <- expand.grid(X, X)
  #   pairs$tag <- pairs[,1] - pairs[,2]
  pairs$multi <- pairs[,1]*pairs[,2]*pop_vector
  return(pairs$multi)
}  
get_cov <- function(K){
  test <- apply(K, MARGIN = 2, FUN = "get_pairs", 
                pop_vector=(nt[[doSpp]]))
  mat_dim <- sqrt(dim(test)[1])
  test <- as.data.frame(test)
  test$tag <- rep(c(1:mat_dim), each=mat_dim)
  cov_str <- matrix(ncol=mat_dim, nrow=mat_dim)
  for(do_i in 1:mat_dim){
    tmp <- subset(test, tag==do_i) #subset out the focal i
    rmtmp <- which(colnames(tmp)=="tag") #get rid of id column
    # Sum over k columns
    cov_str[do_i,] <- (-h[doSpp]^2) * apply(tmp[,-rmtmp], MARGIN = 2, FUN = "sum")
  }
  diag(cov_str) <- 1
  return(cov_str)
}
GenerateMultivariatePoisson<-function(pD, samples, R, lambda){
  normal_mu=rep(0, pD)
  normal = mvrnorm(samples, normal_mu, R)
  pois = normal
  p=pnorm(normal)
  for (s in 1:pD){pois[s]=qpois(p[s], lambda[s])}
  return(pois)
}

make.R.values=function(v,u, #state variables
                       Rpars,rpa,doYear,doSpp){
  f(v,u,Rpars,rpa,doSpp)
}

make.P.values <- function(v,u,muWG,muWS, #state variables
                          Gpars,Spars,doYear,doSpp){  #growth arguments
  S(u,muWS,Spars,doYear,doSpp)*G(v,u,muWG,Gpars,doYear,doSpp) 
}

make.P.matrix <- function(v,muWG,muWS,Gpars,Spars,doYear,doSpp) {
  muWG=expandW(v,v,muWG)
  muWS=expandW(v,v,muWS)
  
  P.matrix=outer(v,v,make.P.values,muWG,muWS,Gpars,Spars,doYear,doSpp)
  return(h[doSpp]*P.matrix)
} 

make.R.matrix=function(v,Rpars,rpa,doYear,doSpp) {
  R.matrix=outer(v,v,make.R.values,Rpars,rpa,doYear,doSpp)
  return(h[doSpp]*R.matrix)
}

# Function to format the W matrix for the outer product
expandW=function(v,u,W){
  if(dim(W)[1]!=length(u)) stop("Check size of W")
  Nspp=dim(W)[2]
  W=as.vector(W)
  W=matrix(W,length(W),ncol=length(v))
  W=as.vector(t(W))
  W=matrix(W,nrow=length(u)*length(v),ncol=Nspp)
  return(W)
}


# Function to calculate size-dependent crowding, assuming no overlap
wrijG=function(r,i,j){
  return(2*pi*integrate(function(z) z*exp(-alphaG[i,j]*(z^2))*Cr[[j]](z-r),r,r+r.U[j])$value+
           pi*Ctot[j]*exp(-alphaG[i,j]*((r+r.U[j])^2))/alphaG[i,j]);   
}
WrijG=Vectorize(wrijG,vectorize.args="r")

wrijS=function(r,i,j){
  return(2*pi*integrate(function(z) z*exp(-alphaS[i,j]*(z^2))*Cr[[j]](z-r),r,r+r.U[j])$value+
           pi*Ctot[j]*exp(-alphaS[i,j]*((r+r.U[j])^2))/alphaS[i,j]);   
}
WrijS=Vectorize(wrijS,vectorize.args="r")


# Function to sum total cover of each species
sumCover=function(v,nt,h,A){
  out=lapply(1:Nspp,function(i,v,nt,h,A) h[i]*sum(nt[[i]]*exp(v[[i]]))/A,v=v,nt=nt,h=h,A=A)
  return(unlist(out))
} 

# Function to sum total density of each species
sumN=function(nt,h){
  out=lapply(1:Nspp,function(i,nt,h) h[i]*sum(nt[[i]]),nt=nt,h=h)
  return(unlist(out))
}

# Function to calculate size variance of each species
varN=function(v,nt,h,Xbar,N){
  out=lapply(1:Nspp,function(i,v,nt,h,Xbar,N) h[i]*sum((exp(v[[i]]-Xbar[i])^2)*nt[[i]])/N[i],v=v,nt=nt,h=h,Xbar=Xbar,N=N)
  return(unlist(out))
}  


####
####  Calculate the equilibrium areas.
####
## initial population density vector
new.nt=nt

# set up matrix to record cover
covSave = matrix(NA,tlimit,Nspp)
covSave[1,]=sumCover(v,nt,h,A)

# set up list to store size distributions
sizeSave=list(NULL)
for(i in 1:Nspp){
  sizeSave[[i]]=matrix(NA,length(v[[i]]),(tlimit))
  sizeSave[[i]][,1]=nt[[i]]/sum(nt[[i]])
}

# initial densities 
Nsave=matrix(NA,tlimit,Nspp)
Nsave[1,]=sumN(nt,h)

yrSave=rep(NA,tlimit)

####
####  Run simulation
####
covmat <- list()
for (i in 2:(tlimit)){
  
  #draw from observed year effects
  allYrs=c(1:Nyrs)
  doYear=sample(allYrs,1)
  yrSave[i]=doYear
  
  #get recruits per area
  cover=covSave[i-1,]; N=Nsave[i-1,]
  rpa=get.rpa(Rpars,cover,doYear)
  
  #calculate size-specific crowding
  alphaG=Gpars$alpha 
  alphaS=Spars$alpha 
  
  
  if(NoOverlap.Inter==F){#T: heterospecific genets cannot overlap; F: overlap allowed
    for(ii in 1:Nspp){ 
      # first do all overlap W's
      Xbar=cover*A/N       # multiply by A to get cover back in cm^2
      varX=varN(v,nt,h,Xbar,N) 
      
      muWG = pi*Xbar*N/(A*alphaG[ii,])
      muWS = pi*Xbar*N/(A*alphaS[ii,])
      
      muWG[is.na(muWG)]=0
      muWS[is.na(muWS)]=0
      
      WmatG[[ii]]=matrix(muWG,nrow=length(v[[ii]]),ncol=Nspp,byrow=T)
      WmatS[[ii]]=matrix(muWS,nrow=length(v[[ii]]),ncol=Nspp,byrow=T)
      
      # now do conspecific no overlap W
      Ctot[ii]=h[ii]*sum(expv[[ii]]*nt[[ii]]) 
      Cr[[ii]]=splinefun(b.r[[ii]],h[ii]*c(0,cumsum(expv[[ii]]*nt[[ii]])),method="natural")
      
      WmatG[[ii]][,ii]=WrijG(v.r[[ii]],ii,ii)/A
      WmatS[[ii]][,ii]=WrijS(v.r[[ii]],ii,ii)/A
    }
  }else{
    for(ii in 1:Nspp){
      Ctot[ii]=h[ii]*sum(expv[[ii]]*nt[[ii]]) 
      Cr[[ii]]=splinefun(b.r[[ii]],h[ii]*c(0,cumsum(expv[[ii]]*nt[[ii]])),method="natural") 
    }
    for(jj in 1:Nspp){ 
      
      WfunG=splinefun(size.range,WrijG(size.range,jj,jj))
      WfunS=splinefun(size.range,WrijS(size.range,jj,jj))
      
      for(ii in 1:Nspp) { 
        WmatG[[ii]][,jj]=WfunG(v.r[[ii]])/A 
        WmatS[[ii]][,jj]=WfunS(v.r[[ii]])/A 
      }
    }
    
  } # end NoOverlap if
  
  for(doSpp in 1:Nspp){  
    if(cover[doSpp]>0){    
      # make kernels and project
      P.matrix <- make.P.matrix(v[[doSpp]],WmatG[[doSpp]],WmatS[[doSpp]],Gpars,Spars,doYear,doSpp)  
      R.matrix <- make.R.matrix(v[[doSpp]],Rpars,rpa,doYear,doSpp)  
      newK <- P.matrix+R.matrix
      covmat[[doSpp]] <- get_cov(K=P.matrix)
      pCont <- GenerateMultivariatePoisson(pD = length(nt[[doSpp]]),
                                           samples = 1,
                                           R = covmat[[doSpp]],
                                           lambda = P.matrix%*%nt[[doSpp]])
      rCont <- rpois(length(nt[[doSpp]]),R.matrix%*%nt[[doSpp]])
      new.nt[[doSpp]] <- pCont+rCont
#       sizeSave[[doSpp]][,i]=new.nt[[doSpp]]/sum(new.nt[[doSpp]])  
    }    
  } # next species
  
#   nt=new.nt 
#   covSave[i,]=sumCover(v,nt,h,A)  # store the cover as cm^2/cm^2
#   Nsave[i,]=sumN(nt,h)
#   
#   print(i)
#   flush.console()
#   
#   if(sum(is.na(nt))>0) browser()  
} # next time step



####
####  Read in IBM covariance matrix for comparison
####
ibm_covmat <- readRDS("../results/ibm_covmat.RDS")
for(i in 1:3) ibm_covmat[[i]][which(is.na(ibm_covmat[[i]])==TRUE)] <-0
for(i in 1:3) ibm_covmat[[i]] <- ibm_covmat[[i]][2:(nrow(ibm_covmat[[i]])),2:(nrow(ibm_covmat[[i]]))]

# pdf("../results/ipm_ibm_covmatdiffs.pdf")
par(mfrow=c(1,3))
for(i in 1:3){
  matdiff <- ibm_covmat[[i]] - covmat[[i+1]]
  hist(matdiff, main=sppList[i+1])
}
library(gplots)
for(i in 1:3){
  matdiff <- ibm_covmat[[i]] - covmat[[i+1]]
  heatmap.2(matdiff,Rowv=FALSE, Colv=FALSE, dendrogram="none",main=sppList[i+1],
            col=cm.colors(20), tracecol="#303030", trace="none", 
            notecol="black", notecex=0.8, keysize = 1.5, margins=c(5, 5))
}
# dev.off()


plot(covmat[[2]], ibm_covmat[[1]], ylim=c(-0.2,0.2), xlim=c(-0.2,0.2))
