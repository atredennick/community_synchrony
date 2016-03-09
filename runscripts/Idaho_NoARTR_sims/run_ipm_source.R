############################################################
##  IPM skeleton to be sourced after set up.
##  
##  Authors: Andrew Tredennick, Peter Adler, Chengjin Chu
##  Email:   atredenn@gmail.com
##  Created: 9-11-2015
############################################################


####
####  Simulation length, Matrix size and initial vectors
####
v=v.r=b.r=expv=Cr=WmatG=WmatS=list(Nspp)
h=r.L=r.U=Ctot=numeric(Nspp)
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


## initial population density vector
nt=v
for(i in 1:Nspp) nt[[i]][]=0.001
if(do_site=="Idaho") nt[[1]] <- 0
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

pb <- txtProgressBar(min=2, max=tlimit, char="+", style=3, width=65)
for (i in 2:(tlimit)){
  
  #draw from observed year effects
  # allYrs=c(1:Nyrs)
  # doYear=sample(allYrs,1)
  # doYear set from sourcing script
  doYear=randyrvec[i]
  yrSave[i]=doYear
  
  #get recruits per area
  cover=covSave[i-1,]; N=Nsave[i-1,]
  rpa=get_rpa(Rpars,cover,doYear,A)
  
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
      K.matrix=make.K.matrix(v[[doSpp]],WmatG[[doSpp]],WmatS[[doSpp]],Rpars,rpa,Gpars,Spars,doYear,doSpp)	
      new.nt[[doSpp]]=K.matrix%*%nt[[doSpp]] 
      sizeSave[[doSpp]][,i]=new.nt[[doSpp]]/sum(new.nt[[doSpp]])  
    }    
  } # next species
  
  nt=new.nt 
  covSave[i,]=sumCover(v,nt,h,A)  # store the cover as cm^2/cm^2
  Nsave[i,]=sumN(nt,h)
  
  setTxtProgressBar(pb, i)
  flush.console()
  
  if(sum(is.na(nt))>0) browser()  
} # next time step
