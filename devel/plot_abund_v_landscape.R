##
##  Plots mean abundance as a function of landscape size
##  From IBM simulations
##

library(plyr)
library(reshape2)

rm(list=ls())

totSims <- 50      # number of simulations per site 
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


output_df <- data.frame(experiment=NA, site=NA, expansion=NA, 
                        run=NA, meannum=NA)

for(do_exp in unique(files_df$experiment)){
  exp_files <- files_df[which(files_df$experiment==do_exp),"filename"]
  
  for(do_site in site_names){ # loop over sites
    tmp_files <- exp_files[grep(do_site, exp_files)]
    nsims <- length(tmp_files)
    
    for(i in 1:nsims){ # loop over plot size simulations
      tmp <- readRDS(as.character(tmp_files[i]))
      tmp <- as.data.frame(tmp, row.names = c(1:nrow(tmp)))
      plantnum.columns <- grep("N.", colnames(tmp))
      
      # Get rid of ARTR column for Idaho, for now
      # if(do_site=="Idaho"){cover.columns <- cover.columns[2:length(cover.columns)]}
      
      end.tmp <- subset(tmp, time==totT)
      extinct <- character(nrow(end.tmp))
      extinct[] <- "no" # creates extinct storage vector 
      
      for(jj in 1:nrow(end.tmp)){ # loop over simulations within plot size
        tmpjj <- subset(end.tmp, run==jj)
        endplantnum <- tmpjj[,plantnum.columns]
        
        # Set extinct to "yes" if no runs coexist
        if(nrow(endplantnum)==0){extinct[jj]<-"yes"}
        
        # Nested if/then for coexistence runs
        if(nrow(endplantnum)>0){
          if(length(which(endplantnum==0))>0){extinct[jj]<-"yes"} 
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
          ts.tmp <- tmpsynch[burn.in:totT,plantnum.columns]
          species_names <- unlist(strsplit(colnames(ts.tmp), "[.]"))
          species_names <- species_names[species_names!="N"]
          num_spp <- ncol(ts.tmp)
          ts.tmp$timestep <- c(1:nrow(ts.tmp))
          colids <- which(colnames(tmpsynch) %in% paste0("N.",species_names))
          mean.plantnum <- colMeans(tmpsynch[burn.in:totT,colids])
          tmpD <- data.frame(experiment=do_exp, site=do_site, expansion=i, 
                             run=count, meannum=mean.plantnum)
          output_df <- rbind(output_df, tmpD)
          count <- count+1
        } # end loop over coexistence runs
        
      } # end if/then for coexistence runs
      
    }# end plot size sim loop
    
  }# end site loop
  
}# end experiment loop

output_df <- output_df[2:nrow(output_df),]
output_df$N <- output_df$meannum*(output_df$expansion^2)
avg_plantnum <- ddply(output_df, .(site, expansion), summarise,
                      value = mean(N))
library(ggplot2)
ggplot(avg_plantnum, aes(x=(expansion^2), y=value, color=site))+
  geom_line()+
  geom_point()

avg_plantnum_oversite <- ddply(avg_plantnum, .(expansion), summarise,
                               avg_value = mean(value))