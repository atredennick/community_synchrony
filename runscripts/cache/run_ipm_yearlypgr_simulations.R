############################################################
##  Script to run IPM simulations for each site.
##  
##  Authors: Andrew Tredennick, Peter Adler, Chengjin Chu
##  Email:   atredenn@gmail.com
##  Created: 9-11-2015
############################################################

##  Clear the workspace
rm(list=ls(all=TRUE))


####
####  Load libraries
####
library(communitySynchrony)
library(boot)
library(mvtnorm)
library(msm)
library(statmod)  


####
####  IPM subroutines
####
## Make combined kernel
make.K.values=function(v,u,muWG,muWS, #state variables
                       Rpars,rpa,Gpars,Spars,doYear,doSpp){  #growth arguments
  f(v,u,Rpars,rpa,doSpp)+S(u,muWS,Spars,doYear,doSpp)*G(v,u,muWG,Gpars,doYear,doSpp) 
}

## Function to make iteration matrix based only on mean crowding
make.K.matrix=function(v,muWG,muWS,Rpars,rpa,Gpars,Spars,doYear,doSpp) {
  muWG=expandW(v,v,muWG)
  muWS=expandW(v,v,muWS)
  K.matrix=outer(v,v,make.K.values,muWG,muWS,Rpars,rpa,Gpars,Spars,doYear,doSpp)
  return(h[doSpp]*K.matrix)
}

## Function to format the W matrix for the outer product
expandW=function(v,u,W){
  if(dim(W)[1]!=length(u)) stop("Check size of W")
  Nspp=dim(W)[2]
  W=as.vector(W)
  W=matrix(W,length(W),ncol=length(v))
  W=as.vector(t(W))
  W=matrix(W,nrow=length(u)*length(v),ncol=Nspp)
  return(W)
}

## Function to calculate size-dependent crowding, assuming no overlap
# Growth
wrijG=function(r,i,j){
  return(2*pi*integrate(function(z) z*exp(-alphaG[i,j]*(z^2))*Cr[[j]](z-r),r,r+r.U[j])$value+
           pi*Ctot[j]*exp(-alphaG[i,j]*((r+r.U[j])^2))/alphaG[i,j]);   
}
WrijG=Vectorize(wrijG,vectorize.args="r")

#Survival
wrijS=function(r,i,j){
  return(2*pi*integrate(function(z) z*exp(-alphaS[i,j]*(z^2))*Cr[[j]](z-r),r,r+r.U[j])$value+
           pi*Ctot[j]*exp(-alphaS[i,j]*((r+r.U[j])^2))/alphaS[i,j]);   
}
WrijS=Vectorize(wrijS,vectorize.args="r")

## Function to sum total cover of each species
sumCover=function(v,nt,h,A){
  out=lapply(1:Nspp,function(i,v,nt,h,A) h[i]*sum(nt[[i]]*exp(v[[i]]))/A,v=v,nt=nt,h=h,A=A)
  return(unlist(out))
} 

## Function to sum total density of each species
sumN=function(nt,h){
  out=lapply(1:Nspp,function(i,nt,h) h[i]*sum(nt[[i]]),nt=nt,h=h)
  return(unlist(out))
}

## Function to calculate size variance of each species
varN=function(v,nt,h,Xbar,N){
  out=lapply(1:Nspp,function(i,v,nt,h,Xbar,N) h[i]*sum(((exp(v[[i]])-Xbar[i])^2)*nt[[i]])/N[i], #the true size 'exp(c[[i]])'
             v=v,nt=nt,h=h,Xbar=Xbar,N=N)
  return(unlist(out))
} 



####
####  Set up simulation parameters
####
Gpars_all <- readRDS("../results/growth_params_list.RDS")
Spars_all <- readRDS("../results/surv_params_list.RDS")
Rpars_all <- readRDS("../results/recruit_parameters.RDS")

randyrlist <- readRDS("../results/randyr_list.RDS")

site_list <- names(Gpars_all)

output_list <- list()

for(do_site in site_list){
  Gpars_tmp <- Gpars_all[[do_site]]
  Spars_tmp <- Spars_all[[do_site]]
  Rpars_tmp <- Rpars_all[[do_site]]
  randyrvec <- randyrlist[[do_site]]
  
  spp_list <- names(Gpars_tmp)
  Nyrs <- nrow(Gpars_tmp[[1]])
  
  #Import and format parameters
  site_path <- paste("../data/", do_site, sep="")
  Gpars <- format_growth_params(do_site = do_site, species_list = spp_list, 
                                Nyrs = Nyrs, Gdata_species = Gpars_tmp)
  Spars <- format_survival_params(do_site = do_site, species_list = spp_list, 
                                  Nyrs = Nyrs, Sdata_species = Spars_tmp)
  Rpars <- format_recruitment_params(do_site = do_site, species_list = spp_list, 
                                     Nyrs = Nyrs, Rdata_species = Rpars_tmp,
                                     path_to_site_data = site_path)
  
  # Set iteration matrix dimensions and max genet sizes by site
  # these are all taken from Chu and Adler 2015 (Ecological Monographs)
  if(do_site=="Arizona"){
    iter_matrix_dims <- c(50,50)
    max_size <- c(170,40)
  }
  if(do_site=="Idaho"){
    iter_matrix_dims <- c(50,75,50,50)
    max_size <- c(3000,202,260,225)
  }
  if(do_site=="Kansas"){
    iter_matrix_dims <- c(75,50,75)
    max_size <- c(1656,823,2056)
  }
  if(do_site=="Montana"){
    iter_matrix_dims <- c(75,50,5,50)
    max_size <- c(2500,130,22,100)
  }
  if(do_site=="NewMexico"){
    iter_matrix_dims <- c(50,50)
    max_size <- c(600,1300)
  }
  
  sim_count <- 1
  site_output <- list()
  
  constant <- FALSE                 # always run with random year effects
  spp_interact <- FALSE             # always run with interspp interactions OFF
  n_spp <- Nspp <- length(spp_list) # all possible forms of Nspp
  A <- 10000                        # area of a 1x1 meter plot, in cm
  maxSize <- max_size               # redo for IPM source script
  bigM <- iter_matrix_dims          # redo for IPM source script
  NoOverlap.Inter <- FALSE          # heterospecifics allowed to overlap
  tlimit <- 2500                    # number of years to simulate
  size_dir <- "../results/stable_size_dists/" # path for stable size dists
  
  ## Run the IPM from a source script
  source("ipm_intrinsic_pgr_source.R") 
  
  ##  Save site output to big list
  output_list[[do_site]] <- maxR[501:tlimit-1,]
  
} # end site loop

## Save the output
saveRDS(output_list, "../results/ipm_yearly_pgr.RDS")

