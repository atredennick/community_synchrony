## Script to plot IBM results

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 9-22-2015

# Clear workspace 
rm(list=ls(all=TRUE))

totSims <- 20      # number of simulations per site (1 here since using large landscape)
totT <- 100      # time steps of simulation
burn.in <- 25    # time steps to discard before calculating cover values
site_colors <- c("grey45", "steelblue", "slateblue4", "darkorange", "purple")

####
####  Load libraries -----------------------------------------------------------
####
library(ggplot2)
library(reshape2)
library(plyr)
library(synchrony)
library(gridExtra)

### Function for inverse of %in%
"%w/o%" <- function(x, y) x[!x %in% y] # x without y


####
####  Read in data and setup up file lists -------------------------------------
####
setwd("../results/ibm_sims/")
all_files <- list.files()
constinterfiles <- all_files[grep("constinter", all_files)]
constnointerfiles <- all_files[grep("constnointer", all_files)]
fluctinterfiles <- all_files[grep("fluctinter", all_files)]
fluctnointerfiles <- all_files[grep("fluctnointer", all_files)]



####
####  Fluctuating environment runs ---------------------------------------------
####
stringlist <- strsplit(fluctinterfiles, split = "_")
site_names <- unique(sapply(stringlist, "[[", 2))
rm(stringlist)


output_df <- data.frame(site=NA, expansion=NA, run=NA, synchrony=NA)
for(do_site in site_names){ # loop over sites
  tmp_files <- fluctinterfiles[grep(do_site, fluctinterfiles)]
  nsims <- length(tmp_files)
  
  for(i in 1:nsims){ # loop over plot size simulations
    tmp <- readRDS(tmp_files[i])
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
        
        synch.run[count] <- as.numeric(community.sync(obs_gr)[1])
        count <- count+1
      } # end loop over coexistence runs
      
      # Create tmp dataframe for rbinding to storage df
      tmpD <- data.frame(synchrony=synch.run) # synchrony of coexistence runs
      tmpD$site <- do_site # current site
      tmpD$expansion <- i # current plot size
      tmpD$run <- c(1:length(runs)) # index for simulation run
      output_df <- rbind(output_df, tmpD[,c("site", "expansion", "run", "synchrony")])
    } # end if/then for coexistence runs
   
  }# end plot size sim loop
  
}# end site loop

# Remove NA first row
output_df <- output_df[2:nrow(output_df),]


### Make plot
fluct.inter.plot <- ggplot(output_df, aes(x=expansion, y=synchrony, color=site))+
  geom_point()+
  stat_smooth(se=FALSE, method="lm", size=1)+
  scale_color_manual(values = site_colors)+
  scale_y_continuous(limits=c(0,1))+
  xlab(expression(paste("Simulated landscape size (", m^2, ")")))+
  ylab("Species synchrony")+
  theme_bw()+
  guides(color=FALSE)+
  ggtitle("A) Fluctuating environment & \ninterspecific interactions")

# output_agg <- ddply(output_df, .(site, expansion), summarise,
#                     avg_synch = mean(synchrony))
# ggplot(output_agg, aes(x=expansion, y=avg_synch, color=site))+
#   geom_line()+
#   geom_point(size=4)




####
####  Constant+inter environment runs ------------------------------------------
####
stringlist <- strsplit(constinterfiles, split = "_")
site_names <- unique(sapply(stringlist, "[[", 2))
rm(stringlist)


output_df <- data.frame(site=NA, expansion=NA, run=NA, synchrony=NA)
for(do_site in site_names){ # loop over sites
  tmp_files <- constinterfiles[grep(do_site, constinterfiles)]
  nsims <- length(tmp_files)
  
  for(i in 1:nsims){ # loop over plot size simulations
    tmp <- readRDS(tmp_files[i])
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
        
        synch.run[count] <- as.numeric(community.sync(obs_gr)[1])
        count <- count+1
      } # end loop over coexistence runs
      
      # Create tmp dataframe for rbinding to storage df
      tmpD <- data.frame(synchrony=synch.run) # synchrony of coexistence runs
      tmpD$site <- do_site # current site
      tmpD$expansion <- i # current plot size
      tmpD$run <- c(1:length(runs)) # index for simulation run
      output_df <- rbind(output_df, tmpD[,c("site", "expansion", "run", "synchrony")])
    } # end if/then for coexistence runs
    
  }# end plot size sim loop
  
}# end site loop

# Remove NA first row
output_df <- output_df[2:nrow(output_df),]


### Make plot
const.inter.plot <- ggplot(output_df, aes(x=expansion, y=synchrony, color=site))+
  geom_point()+
  stat_smooth(se=FALSE, method="lm", size=1)+
  scale_color_manual(values = site_colors)+
  scale_y_continuous(limits=c(0,1))+
  xlab(expression(paste("Simulated landscape size (", m^2, ")")))+
  ylab("Species synchrony")+
  theme_bw()+
  ggtitle("B) Constant environment & \ninterspecific interactions")+
  theme(legend.position = c(0.9, 0.8))



####
####  Fluctuating environment & no interspp comp runs -_------------------------
####
stringlist <- strsplit(fluctnointerfiles, split = "_")
site_names <- unique(sapply(stringlist, "[[", 2))
rm(stringlist)


output_df <- data.frame(site=NA, expansion=NA, run=NA, synchrony=NA)
for(do_site in site_names){ # loop over sites
  tmp_files <- fluctnointerfiles[grep(do_site, fluctnointerfiles)]
  nsims <- length(tmp_files)
  
  for(i in 1:nsims){ # loop over plot size simulations
    tmp <- readRDS(tmp_files[i])
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
        
        synch.run[count] <- as.numeric(community.sync(obs_gr)[1])
        count <- count+1
      } # end loop over coexistence runs
      
      # Create tmp dataframe for rbinding to storage df
      tmpD <- data.frame(synchrony=synch.run) # synchrony of coexistence runs
      tmpD$site <- do_site # current site
      tmpD$expansion <- i # current plot size
      tmpD$run <- c(1:length(runs)) # index for simulation run
      output_df <- rbind(output_df, tmpD[,c("site", "expansion", "run", "synchrony")])
    } # end if/then for coexistence runs
    
  }# end plot size sim loop
  
}# end site loop

# Remove NA first row
output_df <- output_df[2:nrow(output_df),]


### Make plot
fluct.nointer.plot <- ggplot(output_df, aes(x=expansion, y=synchrony, color=site))+
  geom_point()+
  stat_smooth(se=FALSE, method="lm", size=1)+
  scale_color_manual(values = site_colors)+
  scale_y_continuous(limits=c(0,1))+
  xlab(expression(paste("Simulated landscape size (", m^2, ")")))+
  ylab("Species synchrony")+
  theme_bw()+
  guides(color=FALSE)+
  ggtitle("C) Fluctuating environment & \nno interspecific interactions")



####
####  Constant+inter environment runs ------------------------------------------
####
stringlist <- strsplit(constnointerfiles, split = "_")
site_names <- unique(sapply(stringlist, "[[", 2))
rm(stringlist)


output_df <- data.frame(site=NA, expansion=NA, run=NA, synchrony=NA)
for(do_site in site_names){ # loop over sites
  tmp_files <- constnointerfiles[grep(do_site, constnointerfiles)]
  nsims <- length(tmp_files)
  
  for(i in 1:nsims){ # loop over plot size simulations
    tmp <- readRDS(tmp_files[i])
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
        
        synch.run[count] <- as.numeric(community.sync(obs_gr)[1])
        count <- count+1
      } # end loop over coexistence runs
      
      # Create tmp dataframe for rbinding to storage df
      tmpD <- data.frame(synchrony=synch.run) # synchrony of coexistence runs
      tmpD$site <- do_site # current site
      tmpD$expansion <- i # current plot size
      tmpD$run <- c(1:length(runs)) # index for simulation run
      output_df <- rbind(output_df, tmpD[,c("site", "expansion", "run", "synchrony")])
    } # end if/then for coexistence runs
    
  }# end plot size sim loop
  
}# end site loop

# Remove NA first row
output_df <- output_df[2:nrow(output_df),]


### Make plot
const.nointer.plot <- ggplot(output_df, aes(x=expansion, y=synchrony, color=site))+
  geom_point()+
  stat_smooth(se=FALSE, method="lm", size=1)+
  scale_color_manual(values = site_colors)+
  scale_y_continuous(limits=c(0,1))+
  xlab(expression(paste("Simulated landscape size (", m^2, ")")))+
  ylab("Species synchrony")+
  theme_bw()+
  guides(color=FALSE)+
  ggtitle("D) Constant environment & \nno interspecific interactions")




####
####  Arrange plots and save ---------------------------------------------------
####
out.plot <- grid.arrange(fluct.inter.plot, const.inter.plot, 
                         fluct.nointer.plot, const.nointer.plot,
                         nrow=2, ncol=2)


