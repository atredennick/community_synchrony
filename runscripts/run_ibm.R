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


####
####  Set some global variables for the simulation
####
totSims=1
totT=2500     # time steps of simulation
burn.in=50  # time steps to discard before calculating cover values
L=100       # dimension of square quadrat (cm)
expand=2    # 1 = 1x1 m^2, 2 = 2x2m^2, etc
init.cover=c(0,1,1,1)   # in % cover
maxSize=c(8000,500,500,500)
minSize=0.25
do_site="Idaho"
Nyrs=22
myCol=c("black","darkgreen","blue","red")
doGroup=NA  # NA for spatial avg., values 1-6 for a specific group
constant=T        

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
output=matrix(NA,0,3+2*Nspp)
colnames(output)=c("run","time","yrParams",paste("Cov.",sppList,sep=""),
                   paste("N.",sppList,sep=""))  
for(iSim in 1:totSims){ 
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
  lastID=max(plants[,5])
  
  # vectors for cover and density
  N=matrix(0,totT,Nspp)
  N=colSums(matrix(plants[,1],dim(plants)[1],Nspp)==matrix(1:Nspp,dim(plants)[1],Nspp,byrow=T))
  A=rep(0,Nspp)
  for(i in 1:Nspp){
    A[i]=sum(plants[plants[,1]==i,2])/(expand^2*L^2)
  }
  new.N=rep(0,Nspp); new.A=rep(0,Nspp)
  tmp=c(iSim,1,0,A,N)
  output=rbind(output,tmp)
  
  # plot initial conditions
#   pdf("ibm_example.pdf",height=6,width=6)
#   par(mgp=c(2,0.5,0),tcl=-0.2)
#   symbols(x = plants[,3], y = plants[,4], circles = sqrt(plants[,2]/pi),fg=myCol[plants[,1]],
#           xlim=c(0,L*expand),ylim=c(0,L*expand),main ="Time=1",xlab="x",ylab="y",inches=F,lwd=2)
#   
  for(tt in 2:(totT)){
    
    # draw year effects
    doYr=sample(1:22,1)
    nextplants=plants
    
    # distance matrix
    distMat=getDist(plants,L,expand)
    
    # recruitment
    newplants=recruit(Rpars,sizes=plants[,2],spp=plants[,1],doYear=doYr,lastID=lastID,L,expand)
    
    for(ss in 1:Nspp){
      if(N[ss]>0){ # make sure spp ss is not extinct
        
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
    
    nextplants=nextplants[nextplants[,2]>0,]    # remove dead plants 
    nextplants=rbind(nextplants,newplants)     # add recruits
    
    # output cover and density
    A[]=0; N[]=0
    tmp=aggregate(nextplants[,2],by=list(nextplants[,1]),FUN=sum)
    A[tmp[,1]]=tmp[,2]/(expand^2*L^2)
    tmp=aggregate(rep(1,dim(nextplants)[1]),by=list(nextplants[,1]),FUN=sum)
    N[tmp[,1]]=tmp[,2]/(expand^2)
    
    lastID=max(nextplants[,5])
    plants=nextplants
    if(is.matrix(plants)==F) plants=matrix(plants,nrow=1,ncol=length(plants))
    output=rbind(output,c(iSim,tt,doYr,A,N)) 
    
    # save xy coordinates for spatial analysis
    if(tt>burn.in & tt%%10==0){
      if(sum(plants[,1]==1)>4) {
        tmp=cbind(rep(iSim,dim(plants)[1]),rep(tt,dim(plants)[1]),plants)
        outxy=rbind(outxy,tmp)
      }
    }
    print(tt)
  } # next tt
} # next iSim

matplot(output[,2], output[,4:7]*100, type="l")

plot(density(output[,5]), ylim=c(0,250), lwd=2)
c=1
for(i in 6:7){
  c <- c+1
  lines(density(output[,i]), col=myCol[c], lwd=2)
}

community.sync(output[501:2500,4:7])
saveRDS(output, "../results/idaho_ibm_demogstoch_only.RDS")
