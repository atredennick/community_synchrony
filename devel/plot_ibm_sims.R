## Script to plot IBM results

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 9-22-2015

# Clear workspace 
rm(list=ls(all=TRUE))

sim_time <- 2500
burnin <- 501

####
####  Load libraries -----------------------------------------------------------
####
library(ggplot2)
library(reshape2)
library(plyr)
library(synchrony)



####
####  Read in simulation results -----------------------------------------------
####
setwd("../results/ibm_sims/")
all_files <- list.files()
stringlist <- strsplit(all_files, split = "_")
site_names <- unique(sapply(stringlist, "[[", 2))
rm(stringlist)

output_list <- list()

for(do_site in site_names){
  tmp_files <- all_files[grep(do_site, all_files)]
  nsims <- length(tmp_files)
  site_list <- list()
  for(i in 1:nsims){
    tmp <- readRDS(tmp_files[i])
    cols_i_want <- grep("Cov", colnames(tmp))
    site_list[[paste0("sim",i)]] <- tmp[ ,cols_i_want]
  }# end sim loop
  output_list[[do_site]] <- site_list
  rm(site_list)
}# end site loop
rm(tmp_files)
rm(tmp)
rm(all_files)



####
####  Plot all time series -----------------------------------------------------
####
png("/Users/atredenn/Desktop/all_ts.png", width = 8, height = 8, units="in", res=200)
par(mfrow=c(5,5))
for(do_site in site_names){
  sims <- 5
  if(do_site=="Montana") sims <-4
  for(i in 1:sims){
    tmp <- as.matrix(output_list[[do_site]][i])[[1]]
    matplot(1:nrow(tmp),tmp*100, type="l", main=paste0(do_site, i))
  }
}
dev.off()


####
####  Calculate synchrony for each simulation ----------------------------------
####
sim_synchs <- matrix(ncol=length(site_names), nrow=5)
counter <- 1
for(do_site in site_names){
  tmp_list <- output_list[[do_site]]
  site_synch <- numeric(length(tmp_list))
  for(i in 1:length(tmp_list)){
    tmp_ts <- tmp_list[[i]]
    tmp_ts <- as.matrix(tmp_ts[,which(tmp_ts[nrow(tmp_ts),]>0)])
    tmp_ts <- tmp_ts[(nrow(tmp_ts)-30):nrow(tmp_ts),]
    site_synch[i] <- as.numeric(community.sync(tmp_ts)[1])
  }
  sim_synchs[1:length(site_synch),counter] <- site_synch
  counter <- counter+1
  rm(site_synch)
}
sim_synchs <- as.data.frame(sim_synchs)
colnames(sim_synchs) <- site_names
sim_synchs$size <- c(1,2,4,6,8)

####
####  Plot the results ---------------------------------------------------------
####
plot_synchs <- melt(sim_synchs, id.vars = "size")
ggplot(plot_synchs, aes(x=size, y=value, color=variable))+
  geom_line()+
  geom_point()




