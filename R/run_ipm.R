#' Run multi-species Integral Projection Model
#' 
#' @author Andrew Tredennick
#' @param A Area of quadrat (cm X cm).
#' @param tlimit Length of simulation.
#' @param burn_in Number of years to discard before calculating things.
#' @param spp_list Character vector with list of species for focal site.
#' @param Nyrs Number of years represented in the yearly regression coefficients.
#' @param constant TRUE/FALSE flag for simulating in a constant environment (no random year effects).
#' @param iter_matrix_dims The size of the big matrix for each species.
#' @param max_size Maximum size, in cm^2, for allowable for each species (in alphabetical order).
#' @param Rpars Recruitment regression parameters for focal site and species.
#' @param Spars Survival regression parameters for focal site and species.
#' @param Gpars Recruitment regression parameters for focal site and species.

run_ipm <- function(A=10000, tlimit=2500, burn.in=500, sppList,
                    Nyrs, constant=FALSE, bigM, maxSize,
                    Rpars, Spars, Gpars){
  library(IPMdoit)
  library(boot)
  library(mvtnorm)
  library(msm)
  library(statmod)

  # Turn off random year effects if constant==TRUE
  if(constant==TRUE){
    Rpars$intcpt.yr=matrix(Rpars$intcpt.mu,Nyrs,Nspp,byrow=T)
    Gpars$intcpt.yr[]=0;Gpars$slope.yr[]=0
    Spars$intcpt.yr[]=0;Spars$slope.yr[]=0    
  }
  
  
  
  
} # end of function