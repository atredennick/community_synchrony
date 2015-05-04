#' Estimate recruitment regression coefficients using MCMC
#' 
#' @author Andrew Tredennick
#' @param dataframe Merged time series dataframe of recruitment area for all species.
#' @param n_adapt Number of iterations for adaptation of MCMC (default = 5000).
#' @param n_update Number of iterations for update phase of MCMC (default = 10000).
#' @param n_samples Number of samples to collect from MCMC after update phase (default = 20000).
#' @param n_thin Number of iterations by which to thin the samples (default = 50).
#' @param n_chains Number of MCMC chains to sample (default = 3).
#' @param sppList Character vector of species code names for the site.
#' @return Dataframe with named regression coefficients.

recruit_mcmc <- function(dataframe, n_adapt=5000, n_update=10000, 
                         n_samples=20000, n_thin=50, n_chains=3,
                         sppList){
  D <- dataframe
  # Calculate mean cover by group and year
  tmpD <- D[,c("quad","year","Group",paste("cov.",sppList,sep=""))]
  tmpD <- aggregate(tmpD[,4:NCOL(tmpD)],by=list("year"=tmpD$year,"Group"=tmpD$Group),FUN=mean)
  names(tmpD)[3:NCOL(tmpD)] <- paste("Gcov.",sppList,sep="")
  D <- merge(D,tmpD,all.x=T)
  
  # Calculate parent cover
  parents1=as.matrix(D[,c(paste("cov.",sppList,sep=""))])/100 ##convert from absolute cover to [1,100] range
  parents2=as.matrix(D[,c(paste("Gcov.",sppList,sep=""))])/100
  
  # Loop through species to get local and group parents by species
  num_species <- ncol(parents)
  tmp_list <- list()
  for(i in 1:num_species){
    tmpL <- which(parents1[,i]==0) # local
    tmpG <- which(parentsw[,i]==0) # group
    tmp <- intersect(tmpL, tmpG)
    tmp_list[[i]] <- tmp
  } # end species loop
  
  
  tmp <- unique(c(tmp1,tmp2))
} # end function
  



