############################################################
##  Script to run IPM simulations for each site.
##  
##  Authors: Andrew Tredennick, Peter Adler, Chengjin Chu
##  Email:   atredenn@gmail.com
##  Created: 9-11-2015
############################################################

##  Clear the workspace
rm(list=ls(all=TRUE))

##  Set the IPM simulation run time and burn in
tlimit <- 2500
burn_in <- 500

##  Set up looping vectors for spp interactions
inter_comp <- c(TRUE, FALSE)
sim_names <- c("ENVINTER", "ENVNOINTER")


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
  
  ##  Set environment time series to be constant across simulations
  # randyrvec <- sample(1:Nyrs, size = tlimit, replace = TRUE)
  randyrvec <- randyrlist[[do_site]]
  
  sim_count <- 1
  site_output <- list()
  for(do_comp in inter_comp){
    constant <- FALSE                 # always run with random year effects
    spp_interact <- do_comp           # T/F for spp interactions
    n_spp <- Nspp <- length(spp_list) # all possible forms of Nspp
    A <- 10000                        # area of a 1x1 meter plot, in cm
    maxSize <- max_size               # redo for IPM source script
    bigM <- iter_matrix_dims          # redo for IPM source script
    NoOverlap.Inter <- FALSE          # heterospecifics allowed to overlap
    
    # Turn off random year effects if constant==TRUE
    if(constant==TRUE){
      Rpars$intcpt.yr=matrix(Rpars$intcpt.mu,Nyrs,Nspp,byrow=T)
      Gpars$intcpt.yr[]=0;Gpars$slope.yr[]=0
      Spars$intcpt.yr[]=0;Spars$slope.yr[]=0    
    } # end constant environment if/then
    
    # Turn off competition if spp_interact==FALSE
    if(spp_interact==FALSE){
      rnbtmp <- Rpars$dd
      rnbtmp[] <- 0
      diag(rnbtmp) <- diag(Rpars$dd)
      Rpars$dd <- rnbtmp
      
      gnbtmp <- Gpars$nb
      gnbtmp[] <- 0
      diag(gnbtmp) <- diag(Gpars$nb)
      Gpars$nb <- gnbtmp
      
      snbtmp <- Spars$nb
      snbtmp[] <- 0
      diag(snbtmp) <- diag(Spars$nb)
      Spars$nb <- snbtmp    
    } # end interspecific competition if/then
    
    ## Run the IPM from a source script
    source("run_ipm_source.R") 
    
    ##  Formate and save output in list
    colnames(covSave) <- spp_list
    site_output[[sim_names[sim_count]]] <- covSave[(burn_in+1):tlimit,]
    
    ##  Advance counter
    sim_count <- sim_count+1
    
    ##  Save stable size distribution IF spp interact == TRUE
    if(spp_interact==TRUE){
      outfile <- "stable_size.csv"
      for(i in 1:Nspp){
        outdir <- "../results/stable_size_dists/"
        if(file.exists(outdir)==FALSE){dir.create(outdir)}
        filename <- paste0(outdir,do_site,"_",spp_list[i],"_",outfile)
        tmp <- rowMeans(sizeSave[[i]][,(burn_in+1):tlimit])
        output <- data.frame(v[[i]],tmp)
        names(output) <- c("size","freq")
        write.table(output,filename,row.names=F,sep=",")
      } # end species loop
    } # end spp_interact IF/THEN for stable size saving
    
  } # end species interaction loop
  
  ##  Save site output to big list
  output_list[[do_site]] <- site_output
  print(paste("done with ",do_site))
  
} # end site loop

## Save the output
saveRDS(output_list, "../results/ipm_comp_nocomp_sims.RDS")

