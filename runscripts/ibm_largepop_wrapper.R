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
totSims <- 20      # number of simulations per site (1 here since using large landscape)
totT <- 100      # time steps of simulation
burn.in <- 25    # time steps to discard before calculating cover values
L <- 100          # dimension of square quadrat (cm)
doGroup <- NA     # NA for spatial avg., values for a specific group
constant <- FALSE # TRUE for constant env.; FALSE for random year effects

## Looping over different landscape sizes
expand_vec <- c(1,2,3,4,5) # 1 = 1x1 m^2, 2 = 2x2m^2, etc


####
####  Read in regression parameters; list for all sites and species ------------
####
Gpars_all <- readRDS("../results/growth_params_list.RDS")
Spars_all <- readRDS("../results/surv_params_list.RDS")
Rpars_all <- readRDS("../results/recruit_parameters.RDS")

site_names <- names(Gpars_all)
n_sites <- length(site_names)

do_site <- "NewMexico"
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
    init.cover <- c(0,1,1,1) # in percent cover
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
  for(expand in expand_vec){
    source("ibm_skeleton.R")
    ####
    ####  Save site output; raw time series
    ####
    saveRDS(output, paste0("../results/ibm_sims/ibm_", do_site, "_expand", expand, ".RDS"))
    print(paste("Done with", do_site, "expansion", expand))
  } # end expansion landscape loop
  
} # end site loop

# library(ggplot2)
# output <- as.data.frame(output, row.names = c(1:nrow(output)))
# ggplot(output)+
#   geom_line(aes(x=time, y=Cov.BOER))+
#   geom_line(aes(x=time, y=Cov.SPFL))+
#   facet_wrap("run")
# 
# ### Get non-extinction runs
# "%w/o%" <- function(x, y) x[!x %in% y] # x without y
# cover.columns <- grep("Cov.", colnames(output))
# extinctions <- unique(output[which(output[,4]==0 | output[,5]==0), "run"])
# keeps <- output$run %w/o% extinctions
# coexist <- subset(output, run %in% keeps)
# 
# ggplot(coexist)+
#   geom_line(aes(x=time, y=Cov.BOER))+
#   geom_line(aes(x=time, y=Cov.SPFL))+
#   facet_wrap("run")
