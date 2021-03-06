##  Individual-based model with no environmental variance
##  Used to test against the demographic approximation of the IPM

######  FROM PETER #########################################
# Individually-based model for a multiple species,
# with density-dependence and explicit space (random
# spatial pattern)

# Gaussian distance function

# This version uses the negative binomial recruitment function

# Initialized from observed size distributions
############################################################

##  Clear the workspace and load libraries
rm(list=ls())
library(communitySynchrony)
library(synchrony)
library(boot)
library(mvtnorm)
library(msm)


make.P.values <- function(v,u,muWG,muWS, #state variables
                          Gpars,Spars,doYear,doSpp){  #growth arguments
  S(u,muWS,Spars,doYear,doSpp)*G(v,u,muWG,Gpars,doYear,doSpp) 
}

make.P.matrix <- function(v,muWG,muWS,Gpars,Spars,doYear,doSpp) {
  muWG=expandW(v,v,muWG)
  muWS=expandW(v,v,muWS)
  
  P.matrix=outer(v,v,make.P.values,muWG,muWS,Gpars,Spars,doYear,doSpp)
  return(h[dospps]*P.matrix)
} 


####
####  Set some global variables for the simulation
####
totSims=100
totT=150     # time steps of simulation
burn.in=50  # time steps to discard before calculating cover values
L=100       # dimension of square quadrat (cm)
expand=2    # 1 = 1x1 m^2, 2 = 2x2m^2, etc
init.cover=c(0,1,1,1)   # in % cover
# maxSize=c(8000,500,500,500)
maxSize=c(3000,202,260,225) 
minSize=0.2
iter_matrix_dims <- c(50,75,50,75)
do_site="Idaho"
Nyrs=22
myCol=c("black","darkgreen","blue","red")
doGroup=NA  # NA for spatial avg., values 1-6 for a specific group
constant=T        

out.nt <- list()
out.nt[[1]] <- array(0,c(totT,iter_matrix_dims[2], totSims))
out.nt[[2]] <- array(0,c(totT,iter_matrix_dims[3], totSims))
out.nt[[3]] <- array(0,c(totT,iter_matrix_dims[4], totSims))



####
####  Read in regression parameters, subset for Idaho
####  and format
####
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


####
####  Vital rate functions
####
## Growth
grow=function(Gpars,doSpp,doYear,sizes,crowding){
  # crowding and nb are vectors of dim Nspp
  logsizes=log(sizes)
  mu=Gpars$intcpt[doSpp]+Gpars$slope[doSpp]*logsizes+Gpars$nb[doSpp,]%*%crowding
  tmp=which(mu<log(minSize)*1.5)         # we will kill vanishingly small plants...below
  mu[tmp]=log(minSize)                   # truncate tiny sizes (they cause problems in sigma2)
  sigma2=Gpars$sigma2.a[doSpp]*exp(Gpars$sigma2.b[doSpp]*mu)
  out=exp(rnorm(length(sizes),mu,sqrt(sigma2)))
  if(sum(is.na(out))>0) browser()
  out[tmp]=0                             # here's the killing
  out[out>maxSize[doSpp]]=maxSize[doSpp] # truncate big plants
  out
}

## Survival
survive=function(Spars,doSpp,doYear,sizes,crowding){
  logsizes=log(sizes)
  mu=Spars$intcpt[doSpp]+Spars$slope[doSpp]*logsizes+Spars$nb[doSpp,]%*%crowding
  out=inv.logit(mu)
  out=rbinom(length(sizes),1,out)
  out
}

## Recruitment
recruit=function(Rpars,sizes,spp,doYear,lastID,L,expand){
  # dd is a matrix of dim Nspp X Nspp
  # sizes and spp are vectors of the same length (=N plants)
  # calculate total areas
  totArea=aggregate(sizes,by=list("spp"=spp),FUN=sum)
  # put in missing zeros
  tmp=data.frame("spp"=1:length(sppList))
  totArea=merge(totArea,tmp,all.y=T)
  totArea[is.na(totArea)]=0
  totArea=totArea[,2]/((L*expand)^2)*100  # scale to % cover   
  # calculate recruits
  lambda=rep(NA,Nspp) # seed production
  for(i in 1:Nspp){
    lambda[i]=totArea[i]*exp(Rpars$intcpt.mu[i]+totArea%*%Rpars$dd[i,])
  }
  # number of draws from distribution depends on size of landscape
  NN=rnbinom(length(lambda)*expand^2,mu=lambda,size=Rpars$theta)  
  NN=rowSums(matrix(NN,length(lambda),expand^2))      
  x=y=spp=id=size=NULL
  for(i in 1:Nspp){
    if(NN[i]>0){
      #get recruit sizes 
      size=c(size,exp(rnorm(NN[i],Rpars$sizeMean[i],sqrt(Rpars$sizeVar[i]))))
      if(sum(is.na(size))>0) stop("Check recruit sizes")
      #assign random coordinates
      x=c(x,expand*L*runif(NN[i])); y=c(y,expand*L*runif(NN[i]))
      spp=c(spp,rep(i,NN[i]))
      #assign genet ID's
      if(length(id)>0) lastID=max(id)
      id=c(id,(1:NN[i]+lastID))
    }      
  } # next i
  # output
  size[size<minSize]=minSize
  out=cbind(spp,size,x,y,id)
  return(out)
}

####
####  Initialize some things for simulations
# get observed size distribution for initialization
size.init=list()
for(i in 1:Nspp){
  infile=paste("../data/",do_site,"/",sppList[i],"/",sppList[i],"_genet_xy.csv",sep="")
  tmp=read.csv(infile)
  size.init[[i]]=tmp$area  
}



####
####  Define other functions
####
# Crowding function, assumes toroidal landscape
getDist=function(plants,L,expand){
  xdiff=abs(outer(plants[,3],plants[,3],FUN="-"))
  tmp=which(xdiff>((L*expand)/2))
  xdiff[tmp]=(L*expand)-xdiff[tmp]
  ydiff=abs(outer(plants[,4],plants[,4],FUN="-"))
  tmp=which(ydiff>((L*expand)/2))
  ydiff[tmp]=(L*expand)-ydiff[tmp]
  distMat=sqrt(xdiff^2+ydiff^2) 
  distMat[distMat==0]=NA
  return(distMat)
}

getCrowding=function(plants,alpha,distMat){
  if(dim(plants)[1]>1){  
    distMat=exp(-1*alpha[plants[,1]]*distMat^2)
    sizeMat=matrix(plants[,2],dim(plants)[1],dim(plants)[1])
    distSize=distMat*sizeMat
    out=sapply(1:4,function(i,distSize){ 
      colSums(matrix(distSize[plants[,1]==i,],sum(plants[,1]==i),NCOL(distSize)),na.rm=T)},
      distSize=distSize)
    out=t(out)
  }else{
    out=rep(0,Nspp)
  }
  out
}



####
#### Main simulation loop
####
outxy=matrix(NA,0,7)
colnames(outxy)=c("run","t","spp","size","x","y","id")
# output=matrix(NA,0,3+2*Nspp)
# colnames(output)=c("run","time","yrParams",paste("Cov.",sppList,sep=""),
#                    paste("N.",sppList,sep=""))  

# INITIALIZE by drawing from observed sizes until target cover reached
spp=NULL ; size=NULL
for(iSpp in 1:Nspp){
  if(init.cover[iSpp]>0){
    n.init=round(init.cover[iSpp]*L/mean(size.init[[iSpp]])*expand^2)    
    target=init.cover[iSpp]*L*expand^2
    lower=target-0.1*target
    upper=target+0.1*target
    success=F
    while(success==F){
      sizeTry=sample(size.init[[iSpp]],n.init)
      if(sum(sizeTry)>lower & sum(sizeTry)<upper) success=T
    }
    size=c(size,sizeTry)
    spp=c(spp,rep(iSpp,length(sizeTry)))
  }
}
x=runif(length(size),0,L*expand) ; y=runif(length(size),0,L*expand)
id=rep(1:length(size))
plants=cbind(spp,size,x,y,id)
plantsbig=do.call(rbind, replicate(10, cbind(spp,size,x,y,id), simplify = FALSE))
lastID=max(plants[,5])

# vectors for cover and density
N=matrix(0,totT,Nspp)
N=colSums(matrix(plants[,1],dim(plants)[1],Nspp)==matrix(1:Nspp,dim(plants)[1],Nspp,byrow=T))
A=rep(0,Nspp)
for(i in 1:Nspp){
  A[i]=sum(plants[plants[,1]==i,2])/(expand^2*L^2)
}
new.N=rep(0,Nspp); new.A=rep(0,Nspp)
#   tmp=c(iSim,1,0,A,N)
#   output=rbind(output,tmp)

# plot initial conditions
#   pdf("ibm_example.pdf",height=6,width=6)
#   par(mgp=c(2,0.5,0),tcl=-0.2)
#   symbols(x = plants[,3], y = plants[,4], circles = sqrt(plants[,2]/pi),fg=myCol[plants[,1]],
#           xlim=c(0,L*expand),ylim=c(0,L*expand),main ="Time=1",xlab="x",ylab="y",inches=F,lwd=2)
#   
for(tt in 2:(totT)){
  # draw year effects
  doYr=sample(1:22,1)
  
  for(iSim in 1:totSims){
    nextplants=plants
    # distance matrix
    distMat=getDist(plants,L,expand)
    
    # recruitment
    newplants=recruit(Rpars,sizes=plants[,2],spp=plants[,1],doYear=doYr,lastID=lastID,L,expand)
#     newplants=0
    for(ss in 1:Nspp){
      if(N[ss]>0){ # make sure spp ss is not extinct
        ##First small set
        # growth
        W=getCrowding(plants,Gpars$alpha[ss,],distMat)
        newsizes=grow(Gpars,doSpp=ss,doYear=doYr,sizes=plants[,2],crowding=W)
        if(sum(newsizes==Inf)>0) browser()
        if(is.na(sum(newsizes))) browser()
        # survival
        # uses same W as growth
        live=survive(Spars,doSpp=ss,doYear=doYr,sizes=plants[,2],crowding=W)
        # combine growth and survival
        tmp=which(plants[,1]==ss)  # only alter plants of focal spp        
        nextplants[tmp,2]=newsizes[tmp]*live[tmp]   #update with G and S
        
        
      } # end if no plants
    } # next ss 
    
    if(tt<burn.in){
      nextplants=nextplants[nextplants[,2]>0,]    # remove dead plants 
      nextplants=rbind(nextplants,newplants)     # add recruits
    } 
    if(tt>=burn.in) nextplants=nextplants[nextplants[,2]>0,]    # remove dead plants 
    
#     
    
    # output cover and density
    A[]=0; N[]=0
    tmp=aggregate(nextplants[,2],by=list(nextplants[,1]),FUN=sum)
    A[tmp[,1]]=tmp[,2]/(expand^2*L^2)
    tmp=aggregate(rep(1,dim(nextplants)[1]),by=list(nextplants[,1]),FUN=sum)
    N[tmp[,1]]=tmp[,2]/(expand^2)
    
    lastID=max(nextplants[,5])
    
    # Calculate population vector a la IPM nt
    spps <- c(2,3,4)
    for(dospps in spps){
#       if(length(nextplants[which(nextplants[,1]==dospps),])>1)
      tmp.plants <- nextplants[which(nextplants[,1]==dospps),]
      if(is.null(dim(tmp.plants))) {
        print("NULL ONE!!!!!!!")
        tmpmat <- as.matrix(tmp.plants)
        tmp.plants <- data.frame(spp=tmpmat[1,1],size=tmpmat[2,1],x=tmpmat[3,1],y=tmpmat[4,1],id=tmpmat[5,1])
      }
      Low=log(0.2)
      Up=log(maxSize[dospps])*1.1     
      # boundary points b and mesh points y. Note: b chops up the size interval (L-U) into bigM-equal-sized portions.
      bins <- Low+c(0:iter_matrix_dims[dospps])*(Up-Low)/iter_matrix_dims[dospps]
      #   bins <- seq(from = minSize, to = maxSize[dospps], length.out = iter_matrix_dims[dospps])
      v <- 0.5*(bins[1:iter_matrix_dims[dospps]]+bins[2:(iter_matrix_dims[dospps]+1)])
      h <- bins[2]-bins[1]  
      all.bins <- c(1:iter_matrix_dims[dospps])
      tmp.cut <- cut(log(tmp.plants[,2]), breaks = v, labels = FALSE)
      out.nt[[dospps-1]][tt,which(all.bins%in%tmp.cut==TRUE),iSim] <- table(tmp.cut)
    }
    print(c(iSim,tt))
  } # next iSim
  if(tt<burn.in) plants=nextplants
  if(tt>=burn.in) plants=plants
  
  if(is.matrix(plants)==F) plants=matrix(plants,nrow=1,ncol=length(plants))
  #     output=rbind(output,c(iSim,tt,doYr,A,N)) 
  
  # save xy coordinates for spatial analysis
#   if(tt>burn.in & tt%%10==0){
#     if(sum(plants[,1]==1)>4) {
#       tmp=cbind(rep(iSim,dim(plants)[1]),rep(tt,dim(plants)[1]),plants)
#       outxy=rbind(outxy,tmp)
#     }
#   }
}#next tt
  

####
####  Calculate covariance of population vector
####
#Save out.nt for IPM estimates
# saveRDS(out.nt, "../results/nt_popvec_ibm.RDS")
# Sum over sims for nt vectors
library(plyr)
library(reshape2)
ntdf <- melt(out.nt)
colnames(ntdf) <- c("time", "bins", "sims", "genets", "species")

ntagg <- ddply(ntdf, .(time, bins, species), summarise,
               n = sum(genets))
hist(ntagg$n)

out.cov <- list()
for(i in 1:3){
  tmpdf <- subset(ntagg, species==i & time>50)
  cdf <- dcast(tmpdf, time~bins, value.var = "n")
  out <- cor(cdf) #covariance of entire matrix
  diag(out) <- 1
  out.cov[[i]] <- out
}
saveRDS(out.cov, "../results/ibm_covmat.RDS")

spps <- c(2,3,4)
save.nt <- list()
save.nt[[1]] <- numeric(iter_matrix_dims[2])
save.nt[[2]] <- numeric(iter_matrix_dims[3])
save.nt[[3]] <- numeric(iter_matrix_dims[4])
for(i in 1:3) save.nt[[i]][] <- 0
for(dospps in spps){
  tmp.plants <- plants[which(plants[,1]==dospps),]
  Low=log(0.2)
  Up=log(maxSize[dospps])*1.1     
  # boundary points b and mesh points y. Note: b chops up the size interval (L-U) into bigM-equal-sized portions.
  bins <- Low+c(0:iter_matrix_dims[dospps])*(Up-Low)/iter_matrix_dims[dospps]
  #   bins <- seq(from = minSize, to = maxSize[dospps], length.out = iter_matrix_dims[dospps])
  v <- 0.5*(bins[1:iter_matrix_dims[dospps]]+bins[2:(iter_matrix_dims[dospps]+1)])
  h <- bins[2]-bins[1]  
  all.bins <- c(1:iter_matrix_dims[dospps])
  tmp.cut <- cut(log(tmp.plants[,2]), breaks = v, labels = FALSE)
  save.nt[[dospps-1]][which(all.bins%in%tmp.cut==TRUE)] <- table(tmp.cut)
}
saveRDS(save.nt, "../results/nt_popvecSMALL_ibm.RDS")


##  Get BIG nt vector
out.ntbig <- list()
for(i in 1:3){
  tmpdf <- subset(ntagg, species==i & time==50)
  cdf <- dcast(tmpdf, time~bins, value.var = "n")
  out.ntbig[[i]] <- cdf
}
saveRDS(out.ntbig, "../results/nt_popvecBIG_ibm.RDS")

