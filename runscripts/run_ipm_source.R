NoOverlap_Inter=FALSE
  library(IPMdoit)
  library(boot)
  library(mvtnorm)
  library(msm)
  library(statmod)
  A=A
  tlimit=tlimit
  burn_in=burn_in
  spp_list=spp_list
  Nyrs=Nyrs; constant=do_env_const
  iter_matrix_dims=iter_matrix_dims; max_size=max_size
  Rpars=Rpars; Spars=Spars; Gpars=Gpars
  demographic_stochasticity=do_demo_stoch
  spp_interact = spp_interact

  # Turn off random year effects if constant==TRUE
  if(constant==TRUE){
    Rpars$intcpt.yr=matrix(Rpars$intcpt.mu,Nyrs,Nspp,byrow=T)
    Gpars$intcpt.yr[]=0;Gpars$slope.yr[]=0
    Spars$intcpt.yr[]=0;Spars$slope.yr[]=0    
  }

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
  }
  
  # Build initial vectors and matrices
  Nspp <- n_spp <- length(spp_list)
  inits <- make_inits_ipm(n_spp = n_spp, 
                          iter_matrix_dims = iter_matrix_dims, 
                          max_size = max_size)
  
  ##
  ## Run simulation -----------------------------------------------------------
  ##
  # initial population density vector
  nt <- inits$v
  
  # loop through species really quick to set low initial densities
  for(i in 1:n_spp) nt[[i]][]=0.001
  if(do_site=="Idaho")
    nt[[1]][] <- 0
  new.nt <- nt #set initial density vector to be fed into IPM
  
  # set up matrix to record cover
  covSave <- matrix(NA,tlimit,n_spp)
  covSave[1,] <- sum_cover(inits$v,nt,inits$h,A=A)
  
  # set up list to store size distributions
  sizeSave <- list(NULL)
  for(i in 1:n_spp){
    sizeSave[[i]] <- matrix(NA,length(inits$v[[i]]),(tlimit))
    sizeSave[[i]][,1] <- nt[[i]]/sum(nt[[i]])
  }
  
  # initial densities 
  Nsave <- matrix(NA,tlimit,n_spp)
  Nsave[1,] <- sum_N(nt,inits$h)
  
  yrSave <- rep(NA,tlimit)
  
  # Loop through simulation times and iterate population
  pb <- txtProgressBar(min=2, max=tlimit, char="+", style=3, width=65)
  for (t in 2:(tlimit)){
    #draw from observed year effects
    allYrs <- c(1:Nyrs)
    doYear <- sample(allYrs,1)
    yrSave[t] <- doYear
    
    #get recruits per area
    cover <- covSave[t-1,]
    N <- Nsave[t-1,]
    recs_per_area <- get_rpa(Rpars,cover,doYear,A)
    
    #calculate size-specific crowding
    alphaG <- Gpars$alpha 
    alphaS <- Spars$alpha 
    if(NoOverlap_Inter==F){#T: heterospecific genets cannot overlap; F: overlap allowed
      crowd_list <- crowd_overlap(A, N, inits$vt, inits$h, alphaG, alphaS, 
                                  inits$WmatG, inits$WmatS, n_spp, inits$Ctot,
                                  inits$Cr, inits$b.r, inits$expv, inits$r.U, 
                                  inits$v.r, inits$v, cover=cover, nt=nt)
    }else{
      crowd_list <- crowd_no_overlap(A, inits$vt, inits$h, alphaG, alphaS, inits$WmatG, inits$WmatS,
                                     n_spp, inits$Ctot, inits$Cr, inits$b.r, inits$expv, inits$r.U, inits$v.r,
                                     inits$size_range)
    } # end NoOverlap if
    
    for(doSpp in 1:n_spp){  
      if(cover[doSpp]>0){    
        # Make kernels and project
        K_matrix=make_K_matrix(inits$v[[doSpp]],crowd_list$WmatG[[doSpp]],
                               crowd_list$WmatS[[doSpp]],
                               Rpars,recs_per_area,Gpars,Spars,
                               doYear,doSpp,inits$h,
                               demo_stoch=demographic_stochasticity) 
        
        new.nt[[doSpp]]=K_matrix%*%nt[[doSpp]] 
        sizeSave[[doSpp]][,t]=new.nt[[doSpp]]/sum(new.nt[[doSpp]])  
      }    
    } # next species
    
    nt=new.nt 
    covSave[t,]=sum_cover(inits$v,nt,inits$h,A)  # store the cover as cm^2/cm^2
    Nsave[t,]=sum_N(nt,inits$h)
    
    setTxtProgressBar(pb, t)
    flush.console()
    if(sum(is.na(nt))>0) browser()  
  } # next time step
