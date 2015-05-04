#' Estimate recruitment regression coefficients using MCMC
#' 
#' @author Andrew Tredennick
#' @param dataframe Merged time series dataframe of recruitment area for all species.
#' @param n_adapt Number of iterations for adaptation of MCMC (default = 5000).
#' @param n_update Number of iterations for update phase of MCMC (default = 10000).
#' @param n_samples Number of samples to collect from MCMC after update phase (default = 20000).
#' @param n_thin Number of iterations by which to thin the samples (default = 50).
#' @param sppList Character vector of species code names for the site.
#' @return Matrix of statistical results by fitted parameter.

recruit_mcmc <- function(dataframe, n_adapt=5000, n_update=10000, 
                         n_samples=20000, n_thin=50, sppList){
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
  num_species <- length(sppList)
  tmp_list <- list()
  for(i in 1:num_species){
    tmpL <- which(parents1[,i]==0) # local
    tmpG <- which(parents2[,i]==0) # group
    tmp <- intersect(tmpL, tmpG)
    tmp_list[[i]] <- tmp
  } # end species loop
  if(num_species == 2)
    tmp <- unique(c(tmp_list[[1]],tmp_list[[2]]))
  if(num_species == 3)
    tmp <- unique(c(tmp_list[[1]],tmp_list[[2]],tmp_list[[3]]))
  if(num_species == 4)
    tmp <- unique(c(tmp_list[[1]],tmp_list[[2]],tmp_list[[3]],tmp_list[[4]]))
  
  if(length(tmp)>0){
    parents1 <- parents1[-tmp,] ##remove them
    parents2 <- parents2[-tmp,] ##remove them
    y <- as.matrix(D[,c(paste("R.",sppList,sep=""))])[-tmp,] ##remove them  
    year <- as.numeric(as.factor(D$year))[-tmp] ##remove them
    Nyrs <- length(unique(D$year))
    N <- dim(D)[1]-length(tmp) ##reduce
    Nspp <- length(sppList)
    Group <- as.numeric(as.factor(D$Group))[-tmp] ##remove them ##first turn it as FACTOR, then to NUMERIC
    Ngroups <- length(unique(Group))
  } else {
    y <- as.matrix(D[,c(paste("R.",sppList,sep=""))])
    year <- as.numeric(as.factor(D$year))
    Nyrs <- length(unique(D$year))
    N <- dim(D)[1]
    Nspp <- length(sppList)
    Group <- as.numeric(as.factor(D$Group)) ##first turn it as FACTOR, then to NUMERIC
    Ngroups <- length(unique(Group))
  }
  
  # fit as negative binomial with random effects in JAGS
  library(coda)
  library(rjags)
  dataJ = list(N = N, y=y, parents1=parents1, parents2=parents2, 
               year=year, Nyrs=Nyrs, Nspp=Nspp, Ngroups=Ngroups, Group=Group)
  inits <- NULL
  inits[[1]] <- list(intcpt.yr=matrix(1,Nyrs,Nspp),intcpt.mu=rep(1,Nspp),
                     intcpt.tau=rep(1,Nspp),
                     intcpt.gr=matrix(1,Ngroups,Nspp),g.tau=rep(1,Nspp),
                     dd=matrix(-1,Nspp,Nspp),theta=rep(1,Nspp)) 
  inits[[2]] <- list(intcpt.yr=matrix(0,Nyrs,Nspp),intcpt.mu=rep(0,Nspp),
                     intcpt.tau=rep(10,Nspp),
                     intcpt.gr=matrix(0,Ngroups,Nspp),g.tau=rep(0.1,Nspp),
                     dd=matrix(-0.5,Nspp,Nspp),theta=rep(2,Nspp))
  inits[[3]] <- list(intcpt.yr=matrix(0.5,Nyrs,Nspp),intcpt.mu=rep(0.5,Nspp),
                     intcpt.tau=rep(5,Nspp),
                     intcpt.gr=matrix(0.5,Ngroups,Nspp),g.tau=rep(0.2,Nspp),
                     dd=matrix(-0.3,Nspp,Nspp),theta=rep(4,Nspp)) 
  
  params <- c("intcpt.yr","intcpt.mu","intcpt.tau","intcpt.gr",
              "g.tau","dd","theta","u","lambda") 
  
  modelFile <- "../runscripts/recruitJAGS.R"
  
  n.Adapt <- n_adapt
  n.Up <- n_update
  n.Samp <- n_samples
  n.Thin <- n_thin
  
  jm <- jags.model(modelFile, data=dataJ, n.chains=length(inits),
                   inits = inits, n.adapt = n.Adapt)
  update(jm, n.iter=n.Up)
  out <- coda.samples(jm, variable.names=params, n.iter=n.Samp, n.thin=n.Thin)

  zmStat <- summary(out)$stat
  return(zmStat)
} # end function
  



