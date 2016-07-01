#' Define growth kernel function for IPM
#' 
#' @author Andrew Tredennick
#' @param v Population vector.
#' @param u Population vector copy, essentially.
#' @param W Crowding matrix.
#' @param Gpars List of growth parameters.
#' @param doYear Current year for getting year level random effects.
#' @param doSpp Numeric scalar for species.

G_yr <- function(v,u,W,Gpars,doYear,doSpp){
  mu <- Gpars$intcpt[doSpp]+Gpars$intcpt.yr[doYear,doSpp]+
    (Gpars$slope[doSpp]+Gpars$slope.yr[doYear,doSpp])*u+
    W%*%(Gpars$nb[doSpp,]+Gpars$nb.yr[doYear,,doSpp])
  sigma2 <- Gpars$sigma2.a[doSpp]*exp(Gpars$sigma2.b[doSpp]*mu)
  out <- dnorm(v,mu,sqrt(sigma2))
  return(out)
}


#' Define survival kernel function for IPM
#' 
#' @author Andrew Tredennick
#' @param u Population vector.
#' @param W Crowding matrix.
#' @param Spars List of survival parameters.
#' @param doYear Current year for getting year level random effects.
#' @param doSpp Numeric scalar for species.

S_yr <- function(u,W,Spars,doYear,doSpp){
  mu <- Spars$intcpt[doSpp]+Spars$intcpt.yr[doYear,doSpp]+
    (Spars$slope[doSpp]+Spars$slope.yr[doYear,doSpp])*u+
    W%*%(Spars$nb[doSpp,]+Spars$nb.yr[doYear,,doSpp])
  return(inv.logit(mu))
}



#' Define number of recruits per area function.
#' 
#' @author Andrew Tredennick
#' @param Rpars List of recruitment regression parameters.
#' @param cover Current time step cover in m^2. Length is num_species.
#' @param doYear Numeric scalar for current simulation year.

get_rpa_yr=function(Rpars,cover,doYear,A){
  # cover is in m^2 per m^2; convert to % scale:
  cover2=cover*100
  # calculate recruits
  Nspp=length(cover)
  mu=rep(NA,Nspp)
  for(i in 1:Nspp){
    mu[i]=cover2[i]*exp(Rpars$intcpt.yr[doYear,i]+sqrt(cover2)%*%Rpars$dd.yr[doYear,,i]) 
  }
  if(sum(is.na(mu))>0) browser() # stop for errors
  rpa=mu/(cover*A)  # convert from number recruits to recruits per cm^2
  return(rpa)
}

#' Define fecundity function
#' 
#' @author Andrew Tredennick
#' @param v Population vector.
#' @param u Population vector copy, essentially.
#' @param Rpars List of recruitment regression parameters.
#' @param rpa Recruits per area, from get.rpa function.
#' @param doSpp Numeric scalar for current species.

f_yr <- function(v,u,Rpars,rpa,doSpp) { 
  nRecruits <- rpa[doSpp]*exp(u)
  #probability of producing a seedling of size v
  tmp <- dnorm(v,Rpars$sizeMean[doSpp],sqrt(Rpars$sizeVar[doSpp]))/(1-pnorm(-1.61,Rpars$sizeMean[doSpp],sqrt(Rpars$sizeVar[doSpp])))
  #number recruits of each size 
  f <- nRecruits*tmp
  return(f)
}   