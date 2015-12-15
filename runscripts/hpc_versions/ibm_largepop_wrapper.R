##  Wrapper for IBM runs to compare to IPM runs
##
##  The IBM is run on a different landscape sizes (plot size) to
##  to test the effect demographic stochasticity.

##  Author:       Andrew Tredennick, Peter Adler
##  Email:        atredenn@gmail.com
##  Date created: 12-15-2015



####
####  Set Simulation ID
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
sim_id <- as.numeric(myargument)
# sim_id=1

####
####  Load libraries -----------------------------------------------------------
####
library(boot)               # for the IBM
library(mvtnorm)            # for the IBM
library(msm)                # for the IBM
source("import_regression_params_fxns")



####
####  Set up global variables and simulation config ----------------------------
####
totSims <- 50      # number of simulations per site (1 here since using large landscape)
totT <- 100      # time steps of simulation
burn.in <- 25    # time steps to discard before calculating cover values
L <- 100          # dimension of square quadrat (cm)
doGroup <- NA     # NA for spatial avg., values for a specific group
constant.vec <- c(FALSE,TRUE) # TRUE for constant env.; FALSE for random year effects
sppinter.vec <- c(FALSE, TRUE) # TRUE for interspp interactions; FALSE for no interspp interactions
filename.flag <- c("fluctnointer", "constnointer", "fluctinter", "constinter")
sites <- c("Arizona", "Idaho", "Kansas", "Montana", "NewMexico")
## Different landscape sizes
expand_vec <- c(1,2,3,4,5) # 1 = 1x1 m^2, 2 = 2x2m^2, etc

##  Create Simulation Grid
sim.grid <- as.data.frame(expand.grid(constant.vec,sppinter.vec,expand_vec,sites))
colnames(sim.grid) <- c("constant", "sppinter", "size", "site")
sim.grid$fileflag <- rep(rep(filename.flag, times=length(expand_vec)), times=length(sites))



####
####  Read in regression parameters; list for all sites and species ------------
####
Gpars_all <- readRDS("../../results/growth_params_list.RDS")
Spars_all <- readRDS("../../results/surv_params_list.RDS")
Rpars_all <- readRDS("../../results/recruit_parameters.RDS")

site_names <- names(Gpars_all)
n_sites <- length(site_names)



####
####  Set Simulation Parameters
####
constant <- sim.grid[sim_id, "constant"]
sppinter <- sim.grid[sim_id, "sppinter"]
filename.flag.current <- sim.grid[sim_id, "fileflag"]
do_site <- sim.grid[sim_id, "site"]
expand <- sim.grid[sim_id, "size"]
Gpars_site <- Gpars_all[[do_site]]
Spars_site <- Spars_all[[do_site]]
Rpars_site <- Rpars_all[[do_site]]

##  Define bookkeeping variables
sppList <- spp_list <- names(Gpars_site)
Nyrs <- nrow(Gpars_site[[1]])
Nspp <- length(sppList)
site_path <- paste("./data/", do_site, sep="")

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

## Turn off competition if sppinter==FALSE
if(sppinter==FALSE){
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


##  Simulation settings by site (sizes from Chu and Adler 2015)
if(do_site == "Arizona"){
  init.cover <- rep(1, times=Nspp) # in percent cover
  maxSize <- c(170, 40)            # in centimeters
  minSize <- 0.25                  # in centimeters 
} # end Arizona do_site

if(do_site == "Idaho"){
  # init.cover <- c(0,1,1,1)       # in percent cover
  init.cover <- rep(1, times=Nspp) # in percent cover
  maxSize <- c(8000,500,500,500)   # in centimeters
  minSize <- 0.25                  # in centimeters 
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
####  Source ibm_skeleton.R and Save Results
####
source("ibm_skeleton.R")
outfile <- paste0("./results/ibm_sims/ibm_", 
                   do_site, "_", 
                   filename.flag.current, 
                   "expand", expand, ".RDS")
saveRDS(output, outfile)








