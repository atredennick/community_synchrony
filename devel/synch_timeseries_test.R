##  Script to investigate synchrony as function of time-series length
##  For IBM


# Clear workspace 
rm(list=ls(all=TRUE))

totSims <- 100   # number of simulations per site 
totTvec <- seq(10,100,10)      # time steps of simulation
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
# setwd("../results/ibm_sims/")
setwd("../../../../../Volumes/A02046115/ibm_synch/results/")
all_files <- list.files()
fluctnointerfiles <- all_files[grep("constnointer", all_files)]
do_file <- fluctnointerfiles[grep("Idaho", fluctnointerfiles)][5]


tmp <- readRDS(as.character(do_file))
tmp <- as.data.frame(tmp, row.names = c(1:nrow(tmp)))
# cover.columns <- grep("Cov.", colnames(tmp))

output_df <- data.frame(tslength=NA, run=NA, abund_synch=NA)

for(totT in totTvec){
  end.tmp <- subset(tmp, time==totT)
  extinct <- character(nrow(end.tmp))
  extinct[] <- "no" # creates extinct storage vector 
  cover.columns <- grep("Cov.", colnames(tmp))
  
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
      cover.columns <- grep("Cov.", colnames(tmpsynch))
      ts.tmp <- tmpsynch[burn.in:totT,cover.columns]
      species_names <- unlist(strsplit(colnames(ts.tmp), "[.]"))
      species_names <- species_names[species_names!="Cov"]
      num_spp <- ncol(ts.tmp)
      ts.tmp$timestep <- c(1:nrow(ts.tmp))
      colids <- which(colnames(tmpsynch) %in% paste0("Cov.",species_names))
      synch.abund <- as.numeric(community.sync(tmpsynch[burn.in:totT,colids])[1])
      # synch.pgr <- as.numeric(community.sync(obs_gr)[1])
      tmpD <- data.frame(tslength=totT, run=count, abund_synch=synch.abund)
      output_df <- rbind(output_df, tmpD)
      count <- count+1
    } # end loop over coexistence runs
  } # end if/then for coexistence runs
}

matplot(tmpsynch[burn.in:totT,colids],type="l")

output_df <- output_df[2:nrow(output_df),]

aggdf <- ddply(output_df, .(tslength), summarise,
               avg_synch = mean(abund_synch))

boxplot(abund_synch~tslength, data=output_df, ylim=c(0,1), las=1,
        xlab="Length of Time Series (years)", ylab="Synchrony of % Cover",
        main="IBM on 5x5 meter plot")




####
####  Relative Abundances
####
sim_rels <- tmpsynch[,grep("Cov", colnames(tmpsynch))]
sim_rels$totCov <- rowSums(sim_rels) 
sim_avg_rels <- apply(sim_rels[,1:4]/sim_rels[,5],2,mean)
sim_avg_rels2 <- colMeans(sim_rels)[1:4]/colMeans(sim_rels[5])

setwd("~/Repos/community_synchrony/devel")
site_names <- c("Idaho")

for(do_site in site_names){
  path2data <- paste0("../data/", do_site,"/")
  spp_list <- list.files(path2data)
  num_spp <- length(spp_list)
  site_data <- data.frame(quad=NA, year=NA, totCover=NA, species=NA)
  
  for(dospp in 1:num_spp){ #loop through species to read in data
    spp_now <- spp_list[dospp]
    quad_file <- paste0(path2data,spp_now,"/quadratCover.csv")
    spp_data <- read.csv(quad_file)
    spp_data$species <- spp_now
    site_data <- rbind(site_data, spp_data)
  } #end species looping for raw data
  
  site_data <- site_data[2:nrow(site_data),] #remove first NA row
  
  ts_freq <- ddply(site_data, .(year, species), summarise,
                   avg_totcov = mean(totCover))
  ts_freq_tmp <- dcast(ts_freq, year~species, value.var = "avg_totcov")
  ts_freq_tmp$total <- rowSums(ts_freq_tmp[,c(2:ncol(ts_freq_tmp))])
  obs_spp_frequencies <- colMeans(ts_freq_tmp[,c(2:(ncol(ts_freq_tmp)-1))] / ts_freq_tmp[,ncol(ts_freq_tmp)])
}

plot(obs_spp_frequencies, sim_avg_rels2, ylim=c(0,1), xlim=c(0,1))
abline(0,1)




####
####  Look at influence of p's on predictions
####
demo_pred <- function(p) {1/sum(p^(1/2))^2}
p <- matrix(c(seq(0,1,length.out = 100),rev(seq(0,1,length.out = 100))), ncol=2)
poss.preds <- numeric(nrow(p))
for(i in 1:nrow(p)){
  poss.preds[i] <- demo_pred(p[i,])
}

par(mfrow=c(2,1))
plot(p[,1],poss.preds,type="l", xlab="Rel. Abund. of Spp1 (1-Spp2)", 
     ylab="Synchrony", main="Demographic Stochasticity Prediction", las=1)

env_pred <- function(p, env_corrs) {sum(prod(p)*env_corrs)}
p <- matrix(c(seq(0,1,length.out = 100),rev(seq(0,1,length.out = 100))), ncol=2)
poss.preds <- numeric(nrow(p))
env_corrs <- matrix(c(0.5,0.9,0.9,0.5),2,2)
for(i in 1:nrow(p)){
  poss.preds[i] <- env_pred(p[i,],env_corrs)
}

plot(p[,1],poss.preds,type="l", xlab="Rel. Abund. of Spp1 (1-Spp2)", 
     ylab="Synchrony", main="Environmental Stochasticity Prediction", las=1)


