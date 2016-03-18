##  Main plot for synchrony results

library(ggplot2)
library(ggthemes)
library(reshape2)
library(plyr)
library(synchrony)
library(communitySynchrony)

site_colors <- c("purple")


##  Read in IPM results ---------
output_list <- readRDS("noARTR_ipm_comp_nocomp_sims.RDS")
mlist <- melt(output_list)
colnames(mlist)[1:3] <- c("year", "species", "cover")
sites <- unique(mlist$L1)
sims <- unique(mlist$L2)
synch_df <- data.frame(site=NA, experiment=NA, bootnum=NA, 
                       pgr_synch=NA, abund_synch=NA)
boots <- 100
num_iters <- 50
for(dosite in sites){
  tmp_data <- subset(mlist, L1==dosite)
  for(dosim in sims){
    tmpsim <- subset(tmp_data, L2==dosim)[c("year", "species", "cover")]
    for(i in 1:boots){
      begin_year <- sample(x = 1:(max(tmpsim$year)-num_iters), 1)
      end_year <- begin_year+num_iters
      tmp <- subset(tmpsim, year %in% begin_year:end_year)
      tmp <- subset(tmp, species != "ARTR")
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

ipm_synch <- ddply(synch_df, .(site, experiment, typesynch), summarise,
                   mean_synch = mean(synch),
                   up_synch = quantile(synch, 0.95),
                   lo_synch = quantile(synch, 0.05))
ipm_synch <- subset(ipm_synch, typesynch=="pgr_synch")

##  Read in IBM results -----
# Run "collate_ibm_sims.R" first
ibm_synch_all <- readRDS("collated_ibm_sims.RDS") 
ibm_synch_agg <- ddply(ibm_synch_all, .(experiment, site, expansion, typesynch), summarise,
                   avg_synch = mean(synch),
                   up_synch = quantile(synch, 0.95),
                   lo_synch = quantile(synch, 0.05))
ibm_synch <- subset(ibm_synch_agg, expansion==5 & typesynch=="Per capita growth rate")

##  Combine results -----
sim_names <- c("All Drivers", "No D.S.", "No E.S.", "No Comp.", "No Comp. + No D.S.", "No Comp. + No E.S.")
sim_names_order <- paste0(c(3,1,5,4,2,6),sim_names)
site_names <- unique(ipm_synch$site)
nsites <- length(site_names)
site_labels <- site_names
site_labels[which(site_labels=="NewMexico")] <- "New Mexico"

# Get vectors of synchrony for each "experiment"
control_synch <- subset(ibm_synch, experiment=="fluctinter")
nods_synch <- subset(ipm_synch, experiment=="ENVINTER")
noes_synch <- subset(ibm_synch, experiment=="constinter")
nocomp_synch <- subset(ibm_synch, experiment=="fluctnointer")
nodsnocomp_synch <- subset(ipm_synch, experiment=="ENVNOINTER")
noesnocomp_synch <- subset(ibm_synch, experiment=="constnointer")
all_experiments <- c(control_synch$avg_synch,
                     nods_synch$mean_synch,
                     noes_synch$avg_synch,
                     nocomp_synch$avg_synch,
                     nodsnocomp_synch$mean_synch, 
                     noesnocomp_synch$avg_synch)
all_ups <- c(control_synch$up_synch,
             nods_synch$up_synch,
             noes_synch$up_synch,
             nocomp_synch$up_synch,
             nodsnocomp_synch$up_synch, 
             noesnocomp_synch$up_synch)
all_downs <- c(control_synch$lo_synch,
               nods_synch$lo_synch,
               noes_synch$lo_synch,
               nocomp_synch$lo_synch,
               nodsnocomp_synch$lo_synch, 
               noesnocomp_synch$lo_synch)

plot_df <- data.frame(site = rep(site_labels, times=length(sim_names)),
                      simulation = rep(sim_names_order, each=nsites),
                      synchrony = all_experiments,
                      upper_synch = all_ups,
                      lower_synch = all_downs)

##  Make the main (all sims) plot -----
ggplot(plot_df, aes(x=simulation, y=synchrony, fill=site, color=site))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=lower_synch, ymax=upper_synch), color="white", size=1, width=0.25)+
  geom_errorbar(aes(ymin=lower_synch, ymax=upper_synch), width=0.25)+
  scale_fill_manual(values=site_colors, labels=site_labels, name="")+
  scale_color_manual(values=site_colors, labels=site_labels, name="")+
  xlab("Simulation Experiment")+
  ylab("Synchrony of Species' Growth Rates")+
  scale_x_discrete(labels=sim_names)+
  scale_y_continuous(limits=c(0,1))+
  theme_few()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=FALSE,color=FALSE)

ggsave("../../docs/components/all_sims_supp_noARTR.png", width = 4, height = 4, units="in", dpi=75)
  
##  Calculate percent differences for Results
# plot_df$control <- rep(plot_df[which(plot_df$simulation=="1All Drivers"),"synchrony"], times=nrow(plot_df)/nsites)
# plot_df$percent_diff <- with(plot_df, abs(synchrony-control)/((synchrony+control)/2)*100)
# write.csv(plot_df, "../results/synchsims_percent_diffs.csv")

##  Make demographic stochasiticty plot for all landscape sizes -----
# ibm_demo_rms <- subset(ibm_synch_agg, typesynch=="Per capita growth rate")
# ibm_demo_rms <- ibm_demo_rms[which(ibm_demo_rms$experiment %in% c("fluctinter", "constnointer")),]
# 
# ggplot(ibm_demo_rms, aes(x=expansion, y=avg_synch, color=site))+
#   geom_line(aes(linetype=experiment))+
#   geom_point(aes(shape=experiment), size=3)+
#   geom_errorbar(aes(ymin=lo_synch, ymax=up_synch), width=0.25)+
#   scale_color_manual(values=site_colors, labels=site_labels, name="")+
#   xlab(expression(paste("Simulated Area (", m^2,")")))+
#   ylab("Synchrony of Species' Growth Rates")+
#   scale_y_continuous(limits=c(0,1))+
#   scale_shape_discrete(name="",labels=c("No E.S + No Comp.", "All Drivers"))+
#   scale_linetype_discrete(name="",labels=c("No E.S + No Comp.", "All Drivers"))+
#   guides(color=FALSE)+
#   theme_few()+
#   theme(legend.position=c(0.3,0.8),
#         legend.background = element_rect(fill = NA))
# 
# ggsave("../../docs/components/ibm_across_landscape_supp_noARTR.png", width = 4, height = 4, units="in", dpi=75)
# 
# 
# 
