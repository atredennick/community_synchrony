##  Calculate environmental and demographic variances for percent cover
##  and growth rates

library(ggplot2)
library(ggthemes)
library(reshape2)
library(plyr)
library(synchrony)
library(communitySynchrony)

### Function for inverse of %in%
"%w/o%" <- function(x, y) x[!x %in% y] # x without y


####
####  SUPPLEMENTARY FIGURE FOR SPECIES' ENVIRONMENTAL VARIANCES
####  PERCENT COVER
####
##  Read in IPM results ---------
output_list <- readRDS("../results/ipm_comp_nocomp_sims.RDS")
mlist <- melt(output_list)
colnames(mlist)[1:3] <- c("year", "species", "cover")
sites <- unique(mlist$L1)
sims <- unique(mlist$L2)
synch_df <- list()

num_iters <- 50
site_species_variances <- list()
png("../docs/components/environmental_variances_cover.png", width = 10, height=3, units="in", res=100)
par(mfrow=c(1,5), las=3)
for(dosite in sites){
  tmp_data <- subset(mlist, L1==dosite)
  for(dosim in sims){
    tmpsim <- subset(tmp_data, L2==dosim)
    tmp_abund_vars <- matrix(ncol=length(unique(tmpsim$species)), nrow=boots)
    for(i in 1:boots){
      begin_year <- sample(x = 1:(max(tmpsim$year)-num_iters), 1)
      end_year <- begin_year+num_iters
      tmp <- subset(tmpsim, year %in% begin_year:end_year)
      tmpcast <- dcast(tmp, year~species, value.var = "cover")
      tmp_abund_vars[i,] <- (apply(tmpcast[2:ncol(tmpcast)], 2, "sd"))^2
    }# end boots loop
    site_species_variances[[dosite]][[dosim]] <- tmp_abund_vars
    colnames(tmp_abund_vars) <- unique(tmpsim$species)
    if(dosim=="ENVNOINTER") { boxplot(tmp_abund_vars, main=dosite,  ylab="environmental variance", outline=F) }
  }# end experiment/sim loop
}# end site loop
dev.off()



####
####  SUPPLEMENTARY FIGURE FOR SPECIES' ENVIRONMENTAL VARIANCES
####  GROWTH RATES
####
##  Read in IPM results ---------
output_list <- readRDS("../results/ipm_comp_nocomp_sims.RDS")
mlist <- melt(output_list)
colnames(mlist)[1:3] <- c("year", "species", "cover")
sites <- unique(mlist$L1)
sims <- unique(mlist$L2)
synch_df <- list()

num_iters <- 50
boots <- 100
site_species_variances <- list()
png("../docs/components/environmental_variances_pgr.png", width = 10, height=3, units="in", res=100)
par(mfrow=c(1,5), las=3)
for(dosite in sites){
  tmp_data <- subset(mlist, L1==dosite)
  for(dosim in sims){
    tmpsim <- subset(tmp_data, L2==dosim)
    tmp_pgr_vars <- matrix(ncol=length(unique(tmpsim$species)), nrow=boots)
    for(i in 1:boots){
      begin_year <- sample(x = 1:(max(tmpsim$year)-num_iters), 1)
      end_year <- begin_year+num_iters
      tmp <- subset(tmpsim, year %in% begin_year:end_year)
      tmpsynch <- get_ipm_synchrony(tmp)
      tmp_grs <- tmpsynch[["growth_rates"]]
      tmpcast <- dcast(tmp_grs, year~species, value.var = "pgr")
      tmp_pgr_vars[i,] <- (apply(tmpcast[2:ncol(tmpcast)], 2, "sd"))^2
    }# end boots loop
    site_species_variances[[dosite]][[dosim]] <- tmp_pgr_vars
    colnames(tmp_pgr_vars) <- unique(tmpsim$species)
    if(dosim=="ENVNOINTER") { boxplot(tmp_pgr_vars, main=dosite,  ylab="environmental variance", outline=F) }
  }# end experiment/sim loop
}# end site loop
dev.off()



####
####  SUPPLEMENTAL FIGURE FOR SPECIES' DEMOGRAPHIC VARIANCES
####  PERCENT COVER
####
totSims <- 100   # number of simulations per site 
totT <- 100      # time steps of simulation
burn.in <- 25    # time steps to discard before calculating cover values



####
####  Read in data and setup up file lists -------------------------------------
####
# setwd("../results/ibm_sims/")
setwd("../../../../../Volumes/A02046115/ibm_synch/results/")
all_files <- list.files()
constinterfiles <- all_files[grep("constinter", all_files)]
constnointerfiles <- all_files[grep("constnointer", all_files)]
fluctinterfiles <- all_files[grep("fluctinter", all_files)]
fluctnointerfiles <- all_files[grep("fluctnointer", all_files)]
exp_vector <- c("constinter", "constnointer", "fluctinter", "fluctnointer")
files_df <- data.frame(experiment = rep(exp_vector, each=length(constinterfiles)),
                       filename = c(constinterfiles,
                                    constnointerfiles,
                                    fluctinterfiles,
                                    fluctnointerfiles))


####
####  Loop through experiment, sites, and runs ---------------------------------
####

stringlist <- strsplit(fluctinterfiles, split = "_")
site_names <- unique(sapply(stringlist, "[[", 2))
rm(stringlist)


png("../../../../../Users/atredenn/Repos/community_synchrony/docs/components/demographic_variances_cover.png", width = 10, height=3, units="in", res=100)
par(mfrow=c(1,5), las=3)
for(do_exp in "constnointer"){
  exp_files <- files_df[which(files_df$experiment==do_exp),"filename"]
  
  for(do_site in site_names){ # loop over sites
    tmp_files <- exp_files[grep(do_site, exp_files)]
    nsims <- length(tmp_files)
    
    for(i in 5){ # loop over plot size simulations
      tmp <- readRDS(as.character(tmp_files[i]))
      tmp <- as.data.frame(tmp, row.names = c(1:nrow(tmp)))
      cover.columns <- grep("Cov.", colnames(tmp))
      
      # Get rid of ARTR column for Idaho, for now
      # if(do_site=="Idaho"){cover.columns <- cover.columns[2:length(cover.columns)]}
      
      end.tmp <- subset(tmp, time==totT)
      extinct <- character(nrow(end.tmp))
      extinct[] <- "no" # creates extinct storage vector 
      
      for(jj in 1:nrow(end.tmp)){ # loop over simulations within plot size
        tmpjj <- subset(end.tmp, run==jj)
        endcover <- tmpjj[,cover.columns]
        
        # Set extinct to "yes" if no runs coexist
        if(nrow(endcover)==0){extinct[jj]<-"yes"}
        
        # Nested if/then for coexistence runs
        if(nrow(endcover)>0){
          if(length(which(endcover==0))>0){extinct[jj]<-"yes"} 
        } # end if/then for extinction flagging
        
      } # end loop for sims within plot size
      
      # Create tmp dataframe for merging with main data
      ext.df <- data.frame(run=c(1:nrow(end.tmp)), extinct=extinct)
      tmp <- merge(tmp, ext.df)
      tmp4synch <- subset(tmp, extinct=="no")
      
      # If the tmp dataframe is empty, no update
      if(nrow(tmp4synch)==0){output_df <- output_df}
      
      # If the dataframe is not empty, get synchrony for coexistence runs
      if(nrow(tmp4synch)>0){
        runs <- unique(tmp4synch$run)
        synch.run <- numeric(length(runs))
        count <- 1 # set counter for indexing
        
        tmp_vars <- matrix(nrow=length(runs), ncol=length(cover.columns))
        for(k in runs){ # loop over coexistence runs
          tmpsynch <- subset(tmp4synch, run==k)
          ts.tmp <- tmpsynch[burn.in:totT,cover.columns]
          species_names <- unlist(strsplit(colnames(ts.tmp), "[.]"))
          species_names <- species_names[species_names!="Cov"]
          num_spp <- ncol(ts.tmp)
          colids <- which(colnames(tmpsynch) %in% paste0("Cov.",species_names))
          tmp_vars[k,] <- apply(tmpsynch[burn.in:totT,colids], 2, "sd")
        } # end loop over coexistence runs
        
      } # end if/then for coexistence runs
      
    }# end plot size sim loop
    
    colnames(tmp_vars) <- species_names
    boxplot(tmp_vars, main=do_site,  ylab="demographic variance", outline=F)
    
  }# end site loop
  
}# end experiment loop

dev.off()




####
####  SUPPLEMENTAL FIGURE FOR SPECIES' DEMOGRAPHIC VARIANCES
####  GROWTH RATES
####
totSims <- 100   # number of simulations per site 
totT <- 100      # time steps of simulation
burn.in <- 25    # time steps to discard before calculating cover values



####
####  Read in data and setup up file lists -------------------------------------
####
# setwd("../results/ibm_sims/")
setwd("../../../../../Volumes/A02046115/ibm_synch/results/")
all_files <- list.files()
constinterfiles <- all_files[grep("constinter", all_files)]
constnointerfiles <- all_files[grep("constnointer", all_files)]
fluctinterfiles <- all_files[grep("fluctinter", all_files)]
fluctnointerfiles <- all_files[grep("fluctnointer", all_files)]
exp_vector <- c("constinter", "constnointer", "fluctinter", "fluctnointer")
files_df <- data.frame(experiment = rep(exp_vector, each=length(constinterfiles)),
                       filename = c(constinterfiles,
                                    constnointerfiles,
                                    fluctinterfiles,
                                    fluctnointerfiles))


####
####  Loop through experiment, sites, and runs ---------------------------------
####

stringlist <- strsplit(fluctinterfiles, split = "_")
site_names <- unique(sapply(stringlist, "[[", 2))
rm(stringlist)


png("../../../../../Users/atredenn/Repos/community_synchrony/docs/components/demographic_variances_pgr.png", width = 10, height=3, units="in", res=100)
par(mfrow=c(1,5), las=3)
for(do_exp in "constnointer"){
  exp_files <- files_df[which(files_df$experiment==do_exp),"filename"]
  
  for(do_site in site_names){ # loop over sites
    tmp_files <- exp_files[grep(do_site, exp_files)]
    nsims <- length(tmp_files)
    
    for(i in 5){ # loop over plot size simulations
      tmp <- readRDS(as.character(tmp_files[i]))
      tmp <- as.data.frame(tmp, row.names = c(1:nrow(tmp)))
      cover.columns <- grep("Cov.", colnames(tmp))
      
      # Get rid of ARTR column for Idaho, for now
      # if(do_site=="Idaho"){cover.columns <- cover.columns[2:length(cover.columns)]}
      
      end.tmp <- subset(tmp, time==totT)
      extinct <- character(nrow(end.tmp))
      extinct[] <- "no" # creates extinct storage vector 
      
      for(jj in 1:nrow(end.tmp)){ # loop over simulations within plot size
        tmpjj <- subset(end.tmp, run==jj)
        endcover <- tmpjj[,cover.columns]
        
        # Set extinct to "yes" if no runs coexist
        if(nrow(endcover)==0){extinct[jj]<-"yes"}
        
        # Nested if/then for coexistence runs
        if(nrow(endcover)>0){
          if(length(which(endcover==0))>0){extinct[jj]<-"yes"} 
        } # end if/then for extinction flagging
        
      } # end loop for sims within plot size
      
      # Create tmp dataframe for merging with main data
      ext.df <- data.frame(run=c(1:nrow(end.tmp)), extinct=extinct)
      tmp <- merge(tmp, ext.df)
      tmp4synch <- subset(tmp, extinct=="no")
      
      # If the tmp dataframe is empty, no update
      if(nrow(tmp4synch)==0){output_df <- output_df}
      
      # If the dataframe is not empty, get synchrony for coexistence runs
      if(nrow(tmp4synch)>0){
        runs <- unique(tmp4synch$run)
        synch.run <- numeric(length(runs))
        count <- 1 # set counter for indexing
        
        tmp_vars <- matrix(nrow=length(runs), ncol=length(cover.columns))
        for(k in runs){ # loop over coexistence runs
          tmpsynch <- subset(tmp4synch, run==k)
          ts.tmp <- tmpsynch[burn.in:totT,cover.columns]
          species_names <- unlist(strsplit(colnames(ts.tmp), "[.]"))
          species_names <- species_names[species_names!="Cov"]
          num_spp <- ncol(ts.tmp)
          ts.tmp$timestep <- c(1:nrow(ts.tmp))
          
          # Transform to per capita growth rates
          lagts <- ts.tmp
          lagts$lagtimestep <- lagts$timestep+1
          colnames(lagts)[1:(num_spp+1)] <- paste0(colnames(lagts)[1:(num_spp+1)],"_t0") 
          mergedts <- merge(ts.tmp, lagts, by.x="timestep", by.y="lagtimestep")
          transitions <- nrow(mergedts)
          obs_gr <- matrix(nrow=transitions, ncol=num_spp)
          name_ids1 <- which(colnames(mergedts) %in% paste0("Cov.",species_names))
          name_ids2 <- which(colnames(mergedts) %in% paste0("Cov.",species_names, "_t0"))
          for(ii in 1:transitions){
            obs_gr[ii,] <- as.numeric(log(mergedts[ii,name_ids1]/mergedts[ii,name_ids2]))
          }
          tmp_vars[k,] <- apply(obs_gr, 2, "sd")
        } # end loop over coexistence runs
        
      } # end if/then for coexistence runs
      
    }# end plot size sim loop
    
    colnames(tmp_vars) <- species_names
    boxplot(tmp_vars, main=do_site,  ylab="demographic variance", outline=F)
    
  }# end site loop
  
}# end experiment loop

dev.off()

