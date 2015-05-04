#' Run multi-species Integral Projection Model
#' 
#' @author Andrew Tredennick
#' @param A Area of quadrat (cm X cm).
#' @param tlimit Length of simulation.
#' @param burn.in Number of years to discard before calculating things.
#' @param sppList Character vector with list of species for focal site.
#' @param Nyrs Number of years represented in the yearly regression coefficients.
#' @param constant TRUE/FALSE flag for simulating in a constant environment (no random year effects).
#' @param bigM The size of the big matrix for each species.
#' @param maxSize Maximum size, in cm^2, for allowable for each species (in alphabetical order).
#' @param Rpars Recruitment regression parameters for focal site and species.
#' @param Spars Survival regression parameters for focal site and species.
#' @param Gpars Recruitment regression parameters for focal site and species.

run_ipm <- function(A=10000, tlimit=2500, burn.in=500, sppList,
                    Nyrs, constant=FALSE, bigM, maxSize,
                    Rpars, Spars, Gpars){
  library(IPMdoit)
} # end of function