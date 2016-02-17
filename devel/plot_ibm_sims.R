## Script to plot IBM results

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 9-22-2015

# Clear workspace 
rm(list=ls(all=TRUE))

totSims <- 50      # number of simulations per site 
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
                        run=NA, pgr_synch=NA, abund_synch=NA)

for(do_exp in unique(files_df$experiment)){
  exp_files <- files_df[which(files_df$experiment==do_exp),"filename"]
  
  for(do_site in site_names){ # loop over sites
    tmp_files <- exp_files[grep(do_site, exp_files)]
    nsims <- length(tmp_files)
    
    for(i in 1:nsims){ # loop over plot size simulations
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
          
          colids <- which(colnames(tmpsynch) %in% paste0("Cov.",species_names))
          synch.abund <- as.numeric(community.sync(tmpsynch[burn.in:totT,colids])[1])
          synch.pgr <- as.numeric(community.sync(obs_gr)[1])
          tmpD <- data.frame(experiment=do_exp, site=do_site, expansion=i, 
                             run=count, pgr_synch=synch.pgr, 
                             abund_synch=synch.abund)
          output_df <- rbind(output_df, tmpD)
          count <- count+1
        } # end loop over coexistence runs
        
      } # end if/then for coexistence runs
      
    }# end plot size sim loop
    
  }# end site loop

}# end experiment loop



# Remove NA first row
output_df <- output_df[2:nrow(output_df),]
colnames(output_df)[5:6] <- c("Per capita growth rate", "Cover (%)")
synch_df <- melt(output_df, id.vars = c("experiment", "site", "expansion", "run"))
colnames(synch_df) <- c("experiment", "site", "expansion", "run", 
                        "typesynch", "synch")

# Check for single datapoints
allcombo <- ddply(synch_df, .(experiment, site, expansion, typesynch), summarise,
                  nums = length(run))
singles <- allcombo[which(allcombo$nums<10),]
singles$remove="yes"
test <- merge(synch_df, singles, by=c("experiment", "site", "expansion", "typesynch"), all = TRUE)
ids.out <- which(test$remove=="yes")
synch_dfplot <- test[-ids.out,]
avg_synchs <- ddply(synch_dfplot, .(experiment, expansion, typesynch), summarise,
                    avg_synch = mean(synch))

mono_poly <- subset(synch_dfplot, experiment=="fluctinter" & expansion==5)
saveRDS(mono_poly, file = "../../../../../Users/atredenn/Repos/community_synchrony/devel/mono_poly_ibm.RDS")


### Make plot and save
ibm_plot <- ggplot(synch_dfplot, aes(x=expansion, y=synch, color=experiment))+
  geom_point(alpha=0.5)+
  stat_smooth(se=FALSE, method="lm", size=0.7)+
  scale_color_manual(values = c("steelblue", "slateblue4", "darkorange", "darkred"),
                     name="Simulation \nType",
                     labels=c("C-INTER", "C-NoINTER",
                              "F-INTER", "F-NoINTER"))+
  scale_y_continuous(limits=c(0,1))+
  xlab(expression(paste("Simulated landscape size (", m^2, ")")))+
  ylab("Community synchrony")+
  facet_grid(typesynch~site)+
  theme_bw()

setwd("~/Repos/community_synchrony/devel")
png("../docs/components/ibm_sims_fig.png", height = 5, 
    width = 10, units = "in", res=150)
print(ibm_plot)
dev.off()

