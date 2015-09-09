##  Wrapper for IBM runs to compare to IPM runs
##
##  The IBM is run on a large enough landscape (plot size) to
##  create populations that are "immune" to demographic stochasticity.
##  We can then compare these simulations to those from the IPM to see
##  if we get similar results. In theory (Ellner IPM book), we should.



# Clear the workspace
rm(list=ls())

####
####  Load libraries -----------------------------------------------------------
####
library(communitySynchrony) # inhouse package
library(synchrony)          # calculates Claire's synchrony metric
library(boot)               # for the IBM
library(mvtnorm)            # for the IBM
library(msm)                # for the IBM



####
####  Set up global variables and simulation config ----------------------------
####
totSims <- 1      # number of simulations per site (1 here since using large landscape)
totT <- 2500      # time steps of simulation
burn.in <- 500    # time steps to discard before calculating cover values
L <- 100          # dimension of square quadrat (cm)
expand <- 1       # 1 = 1x1 m^2, 2 = 2x2m^2, etc
doGroup <- NA     # NA for spatial avg., values for a specific group
constant <- FALSE # TRUE for constant env.; FALSE for random year effects



####
####  Read in regression parameters; list for all sites and species ------------
####
Gpars_all <- readRDS("../results/growth_params_list.RDS")
Spars_all <- readRDS("../results/surv_params_list.RDS")
Rpars_all <- readRDS("../results/recruit_parameters.RDS")

site_names <- names(Gpars_all)
n_sites <- length(site_names)


####
####  Start loop over sites ----------------------------------------------------  
####
for(do_site in site_names){
  ##  Get site-specific regression parameters
  Gpars_site <- Gpars_all[[do_site]]
  Spars_site <- Spars_all[[do_site]]
  Rpars_site <- Rpars_all[[do_site]]
  
  ##  Define bookkeeping variables
  sppList <- spp_list <- names(Gpars_site)
  Nyrs <- nrow(Gpars_site[[1]])
  Nspp <- length(sppList)
  site_path <- paste("../data/", do_site, sep="")
  
  ##  Formate regression parameters for model
  Gpars <- format_growth_params(do_site = do_site, species_list = spp_list, 
                                Nyrs = Nyrs, Gdata_species = Gpars_site)
  Spars <- format_survival_params(do_site = do_site, species_list = spp_list, 
                                  Nyrs = Nyrs, Sdata_species = Spars_site)
  Rpars <- format_recruitment_params(do_site = do_site, species_list = spp_list, 
                                     Nyrs = Nyrs, Rdata_species = Rpars_site,
                                     path_to_site_data = site_path)
  
  ## Turn off random year effects if constant==TRUE
  if(constant==TRUE){
    Rpars$intcpt.yr <- matrix(Rpars$intcpt.mu,Nyrs,Nspp,byrow=T)
    Gpars$intcpt.yr[] <- 0; Gpars$slope.yr[] <- 0
    Spars$intcpt.yr[] <- 0; Spars$slope.yr[] <- 0    
  } # end constant T/F
  
  
  ##  Simulation settings by site (sizes from Chu and Adler 2015)
  if(do_site == "Arizona"){
    init.cover <- rep(1, times=Nspp) # in percent cover
    maxSize <- c(170, 40)            # in centimeters
    minSize <- 0.25                  # in centimeters 
  } # end Arizona do_site
  
  if(do_site == "Idaho"){
    init.cover <- c(0,1,1,1)       # in percent cover
    maxSize <- c(8000,500,500,500) # in centimeters
    minSize <- 0.25                # in centimeters 
  } # end Idaho do_site
  
  if(do_site == "Kansas"){
    init.cover <- rep(1, times=Nspp) # in percent cover
    maxSize <- c(1650, 550, 2056)    # in centimeters
    minSize <- 0.25                  # in centimeters 
  } # end Kansas do_site
  
  if(do_site == "Montana"){
    init.cover <- rep(1, times=Nspp) # in percent cover
    maxSize <- c(2500, 130, 22, 100) # in centimeters
    minSize <- 0.25                  # in centimeters 
  } # end Montana do_site
  
  if(do_site == "NewMexico"){
    init.cover <- rep(1, times=Nspp) # in percent cover
    maxSize <- c(600, 1300)          # in centimeters
    minSize <- 0.25                  # in centimeters 
  } # end New Mexico do_site
  
  
  ####
  ####  Source ibm_skeleton.R
  ####
  source("ibm_skeleton.R")

  ##  Make sure nothing went extinct
  cov_cols <- grep("Cov", colnames(output))
  if(0 %in% output[totT,cov_cols]){
    stop("at least one species went extinct, try a larger landscape")
  } # end extinction error IF/THEN
  
  
  ####
  ####  Save site output; raw time series
  ####
  saveRDS(output[(burn.in+1):totT, ], paste0("../results/ibm_largepop_", do_site, ".RDS"))
  
  print("")
  print(paste("Done with", do_site))
  
} # end site loop


