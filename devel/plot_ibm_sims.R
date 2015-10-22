## Script to plot IBM results

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 9-22-2015

# Clear workspace 
rm(list=ls(all=TRUE))

totSims <- 20      # number of simulations per site (1 here since using large landscape)
totT <- 100      # time steps of simulation
burn.in <- 25    # time steps to discard before calculating cover values

####
####  Load libraries -----------------------------------------------------------
####
library(ggplot2)
library(reshape2)
library(plyr)
library(synchrony)


"%w/o%" <- function(x, y) x[!x %in% y] # x without y


####
####  Read in simulation results -----------------------------------------------
####
setwd("../results/ibm_sims/")
all_files <- list.files()
stringlist <- strsplit(all_files, split = "_")
site_names <- unique(sapply(stringlist, "[[", 2))
rm(stringlist)


output_df <- data.frame(site=NA, expansion=NA, run=NA, synchrony=NA)
for(do_site in site_names){
  tmp_files <- all_files[grep(do_site, all_files)]
  nsims <- length(tmp_files)
  for(i in 1:nsims){
    tmp <- readRDS(tmp_files[i])
    tmp <- as.data.frame(tmp, row.names = c(1:nrow(tmp)))
    cover.columns <- grep("Cov.", colnames(tmp))
    if(do_site=="Idaho"){cover.columns <- cover.columns[2:length(cover.columns)]}
    end.tmp <- subset(tmp, time==totT)
    extinct <- character(nrow(end.tmp))
    extinct[] <- "no"
    for(jj in 1:totSims){
      tmpjj <- subset(end.tmp, run==jj)
      endcover <- tmpjj[,cover.columns]
      if(nrow(endcover)==0){extinct[jj]<-"yes"}
      if(nrow(endcover)>0){
        if(length(which(endcover==0))>0){extinct[jj]<-"yes"} 
      }
    }
    ext.df <- data.frame(run=c(1:totSims), extinct=extinct)
    tmp <- merge(tmp, ext.df)
    tmp4synch <- subset(tmp, extinct=="no")
    if(nrow(tmp4synch)==0){output_df <- output_df}
    
    if(nrow(tmp4synch)>0){
      runs <- unique(tmp4synch$run)
      synch.run <- numeric(length(runs))
      count <- 1
      for(k in runs){
        tmpsynch <- subset(tmp4synch, run==k)
        ts.tmp <- as.matrix(tmpsynch[burn.in:totT,cover.columns])
        synch.run[count] <- as.numeric(community.sync(ts.tmp)[[1]])
        count <- count+1
      }
      tmpD <- data.frame(synchrony=synch.run)
      tmpD$site <- do_site
      tmpD$expansion <- i
      tmpD$run <- c(1:length(runs))
      output_df <- rbind(output_df, tmpD[,c("site", "expansion", "run", "synchrony")])
    }
   
  }# end sim loop
}# end site loop
output_df <- output_df[2:nrow(output_df),]

ggplot(output_df, aes(x=expansion, y=synchrony, color=site))+
  geom_point()

output_agg <- ddply(output_df, .(site, expansion), summarise,
                    avg_synch = mean(synchrony))
ggplot(output_agg, aes(x=expansion, y=avg_synch, color=site))+
  geom_line()+
  geom_point(size=4)




####
####  Constant environment runs
####
all_files <- list.files()
dofiles <- grep("constant", all_files)
all_files <- all_files[dofiles]
stringlist <- strsplit(all_files, split = "_")
site_names <- unique(sapply(stringlist, "[[", 2))
rm(stringlist)


output_df <- data.frame(site=NA, expansion=NA, run=NA, synchrony=NA)
for(do_site in site_names){
  tmp_files <- all_files[grep(do_site, all_files)]
  nsims <- length(tmp_files)
  for(i in 1:nsims){
    tmp <- readRDS(tmp_files[i])
    tmp <- as.data.frame(tmp, row.names = c(1:nrow(tmp)))
    cover.columns <- grep("Cov.", colnames(tmp))
    if(do_site=="Idaho"){cover.columns <- cover.columns[2:length(cover.columns)]}
    end.tmp <- subset(tmp, time==totT)
    extinct <- character(nrow(end.tmp))
    extinct[] <- "no"
    for(jj in 1:totSims){
      tmpjj <- subset(end.tmp, run==jj)
      endcover <- tmpjj[,cover.columns]
      if(nrow(endcover)==0){extinct[jj]<-"yes"}
      if(nrow(endcover)>0){
        if(length(which(endcover==0))>0){extinct[jj]<-"yes"} 
      }
    }
    ext.df <- data.frame(run=c(1:totSims), extinct=extinct)
    tmp <- merge(tmp, ext.df)
    tmp4synch <- subset(tmp, extinct=="no")
    if(nrow(tmp4synch)==0){output_df <- output_df}
    
    if(nrow(tmp4synch)>0){
      runs <- unique(tmp4synch$run)
      synch.run <- numeric(length(runs))
      count <- 1
      for(k in runs){
        tmpsynch <- subset(tmp4synch, run==k)
        ts.tmp <- as.matrix(tmpsynch[burn.in:totT,cover.columns])
        synch.run[count] <- as.numeric(community.sync(ts.tmp)[[1]])
        count <- count+1
      }
      tmpD <- data.frame(synchrony=synch.run)
      tmpD$site <- do_site
      tmpD$expansion <- i
      tmpD$run <- c(1:length(runs))
      output_df <- rbind(output_df, tmpD[,c("site", "expansion", "run", "synchrony")])
    }
    
  }# end sim loop
}# end site loop
output_df <- output_df[2:nrow(output_df),]

ggplot(output_df, aes(x=expansion, y=synchrony, color=site))+
  geom_point()

output_agg <- ddply(output_df, .(site, expansion), summarise,
                    avg_synch = mean(synchrony))
ggplot(output_agg, aes(x=expansion, y=avg_synch, color=site))+
  geom_line()+
  geom_point(size=4)




####
####  Plot all time series -----------------------------------------------------
####
# png("/Users/atredenn/Desktop/all_ts.png", width = 8, height = 8, units="in", res=200)
# par(mfrow=c(5,5))
# for(do_site in site_names){
#   sims <- 5
#   if(do_site=="Montana") sims <-4
#   for(i in 1:sims){
#     tmp <- as.matrix(output_list[[do_site]][i])[[1]]
#     matplot(1:nrow(tmp),tmp*100, type="l", main=paste0(do_site, i))
#   }
# }
# dev.off()
# 
# 
# ####
# ####  Calculate synchrony for each simulation ----------------------------------
# ####
# sim_synchs <- matrix(ncol=length(site_names), nrow=5)
# counter <- 1
# for(do_site in site_names){
#   tmp_list <- output_list[[do_site]]
#   site_synch <- numeric(length(tmp_list))
#   for(i in 1:length(tmp_list)){
#     tmp_ts <- tmp_list[[i]]
#     tmp_ts <- as.matrix(tmp_ts[,which(tmp_ts[nrow(tmp_ts),]>0)])
#     tmp_ts <- tmp_ts[(nrow(tmp_ts)-30):nrow(tmp_ts),]
#     site_synch[i] <- as.numeric(community.sync(tmp_ts)[1])
#   }
#   sim_synchs[1:length(site_synch),counter] <- site_synch
#   counter <- counter+1
#   rm(site_synch)
# }
# sim_synchs <- as.data.frame(sim_synchs)
# colnames(sim_synchs) <- site_names
# sim_synchs$size <- c(1,2,4,6,8)
# 
# ####
# ####  Plot the results ---------------------------------------------------------
# ####
# plot_synchs <- melt(sim_synchs, id.vars = "size")
# ggplot(plot_synchs, aes(x=size, y=value, color=variable))+
#   geom_line()+
#   geom_point()
# 
# 
# 
# 
