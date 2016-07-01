######################################################################################
##  plot_montana_ipm_yearlycomp.R: script to calculate and plot the synchrony       ##
##  of species in the Montana community with and without yearly competition effects ##
######################################################################################

rm(list=ls(all.names = TRUE)) # clear the workspace

library(reshape2)
library(plyr)
library(ggplot2)
library(communitySynchrony)
library(synchrony)
library(ggthemes)



##  Read in IPM results ---------
output_list <- readRDS("../results/ipm_comp_nocomp_sims_yearlycomp_Montana.RDS")
mlist <- melt(output_list)
colnames(mlist)[1:3] <- c("year", "species", "cover")
sites <- unique(mlist$L1)
sims <- unique(mlist$L2)
synch_df <- data.frame(site=NA, experiment=NA, bootnum=NA, 
                       pgr_synch=NA, abund_synch=NA)
boots <- 100
num_iters <- 50
for(dosite in "Montana"){
  tmp_data <- subset(mlist, L1==dosite)
  for(dosim in sims){
    tmpsim <- subset(tmp_data, L2==dosim)[c("year", "species", "cover")]
    for(i in 1:boots){
      begin_year <- sample(x = 1:(max(tmpsim$year)-num_iters), 1)
      end_year <- begin_year+num_iters
      tmp <- subset(tmpsim, year %in% begin_year:end_year)
      tmpsynch <- get_ipm_synchrony(tmp)
      tmp_pgr_synch <- as.numeric(tmpsynch$pgr_synchrony["obs"])
      tmpcast <- dcast(tmp, year~species, value.var = "cover")
      tmp_abund_synch <- as.numeric(community.sync(tmpcast[2:ncol(tmpcast)])[1])
      tmp_df <- data.frame(site=dosite, experiment=dosim, 
                           bootnum=i, pgr_synch=tmp_pgr_synch, 
                           abund_synch=tmp_abund_synch)
      synch_df <- rbind(synch_df, tmp_df)
    }# end boots loop
  }# end experiment/sim loop
}# end site loop

synch_dftmp <- synch_df[2:nrow(synch_df),]
synch_df <- melt(synch_dftmp, id.vars = c("site", "experiment", "bootnum"))
colnames(synch_df) <- c("site", "experiment", "bootnum", "typesynch", "synch")

ipm_synch_all <- ddply(synch_df, .(site, experiment, typesynch), summarise,
                       mean_synch = mean(synch),
                       up_synch = quantile(synch, 0.95),
                       lo_synch = quantile(synch, 0.05))
ipm_synch <- subset(ipm_synch_all, typesynch=="pgr_synch")
ipm_synch$experiment <- c("ANo D.S.", "BNo Comp. + No D.S.")

site_colors <- c("grey45", "steelblue", "slateblue4", "darkorange", "purple")

ggplot(ipm_synch, aes(x=experiment, y=mean_synch))+
  geom_bar(stat="identity", fill="darkorange")+
  geom_errorbar(aes(ymin=lo_synch, ymax=up_synch), width=0.25, color="white", size=1)+
  geom_errorbar(aes(ymin=lo_synch, ymax=up_synch), width=0.25, color="darkorange")+
  xlab("Simulation Experiment")+
  ylab("Synchrony of Species' Growth Rates")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_discrete(labels=c("No D.S.", "No Comp. + No D.S."))+
  theme_few()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=FALSE,color=FALSE)
ggsave("../docs/components/Montana_yearlycomp_ipmsynch.png", width = 3, height = 5, units = "in", dpi = 120)
