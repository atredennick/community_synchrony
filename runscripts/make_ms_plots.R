##  Main plot for synchrony results

library(ggplot2)
library(ggthemes)
library(reshape2)
library(plyr)
library(synchrony)
library(communitySynchrony)

site_colors <- c("grey45", "steelblue", "slateblue4", "darkorange", "purple")



####
####  MAIN TEXT FIGURES 1-2
####
##  Read in IPM results ---------
output_list <- readRDS("../results/ipm_comp_nocomp_sims.RDS")
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

##  Read in IBM results -----
ibm_synch_all <- readRDS("../results/ibm_sims_collated.RDS") 
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
site_labels_order <- paste0(c(2,5,3,4,1), site_labels)

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

plot_df <- data.frame(site = rep(site_labels_order, times=length(sim_names)),
                      simulation = rep(sim_names_order, each=nsites),
                      synchrony = all_experiments,
                      upper_synch = all_ups,
                      lower_synch = all_downs)


site_labels <- site_labels[order(site_labels_order)]
sim_labels <- sim_names[order(sim_names_order)]

# Bring in theoretical predictions and observed values
theory_env <- readRDS("../results/envstoch_predictions.RDS")
theory_env <- theory_env[2:nrow(theory_env), c("site", "growthrate_prediction")]
metrics <- readRDS("../results/calculate_site_metrics.RDS")
theory_demo <- metrics[,c("demo_only_exp_grsynch", "obs_gr_synch", "site")]
theory_demo[which(theory_demo$site=="New Mexico"),"site"] <- "NewMexico"
theoretical_preds_and_obs <- merge(theory_env, theory_demo)
colnames(theoretical_preds_and_obs) <- c("site", "envonly_pred", "demonly_pred", "obs")
theoretical_preds_and_obs$site <- site_labels_order

##  Make the main (all sims) plot -----
pointocols <- rep("grey25",3)
ggplot(data=plot_df)+
  geom_bar(data=plot_df, aes(x=simulation, y=synchrony, fill=site, color=site),
           stat="identity")+
  geom_errorbar(data=plot_df, aes(x=simulation, y=synchrony, fill=site, color=site, 
                             ymin=lower_synch, ymax=upper_synch), 
                color="white", size=1, width=0.25)+
  geom_errorbar(data=plot_df, aes(x=simulation, y=synchrony, fill=site, color=site, 
                             ymin=lower_synch, ymax=upper_synch), width=0.25)+
  geom_point(data=theoretical_preds_and_obs, aes(x=3, y=obs), size=3.5, color="white", shape=15)+
  geom_point(data=theoretical_preds_and_obs, aes(x=3, y=obs), size=3, color=pointocols[1], shape=15)+
  geom_point(data=theoretical_preds_and_obs, aes(x=2, y=envonly_pred), size=3.5, color="white", shape=16)+
  geom_point(data=theoretical_preds_and_obs, aes(x=2, y=envonly_pred), size=3, color=pointocols[2], shape=16)+
  geom_point(data=theoretical_preds_and_obs, aes(x=6, y=demonly_pred), size=3.5, color="white", shape=17)+
  geom_point(data=theoretical_preds_and_obs, aes(x=6, y=demonly_pred), size=3, color=pointocols[3], shape=17)+
  facet_wrap("site", nrow=1)+
  scale_fill_manual(values=site_colors, labels=site_labels, name="")+
  scale_color_manual(values=site_colors, labels=site_labels, name="")+
  xlab("Simulation Experiment")+
  ylab("Synchrony of Species' Growth Rates")+
  scale_x_discrete(labels=sim_labels)+
  scale_y_continuous(limits=c(0,1))+
  theme_few()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=FALSE,color=FALSE)

ggsave("../docs/components/all_sims_results.png", width = 10, height = 4, units="in", dpi=75)
  
##  Calculate percent differences for Results
plot_df$control <- rep(plot_df[which(plot_df$simulation=="3All Drivers"),"synchrony"], times=nrow(plot_df)/nsites)
plot_df$percent_diff <- with(plot_df, abs(synchrony-control)/((synchrony+control)/2)*100)
write.csv(plot_df, "../results/synchsims_percent_diffs.csv")

##  Save the plot data
saveRDS(plot_df, "../results/allsims_plot_data.RDS")

##  Make demographic stochasiticty plot for all landscape sizes -----
ibm_demo_rms <- subset(ibm_synch_agg, typesynch=="Per capita growth rate")
ibm_demo_rms <- ibm_demo_rms[which(ibm_demo_rms$experiment %in% c("fluctnointer", "constnointer")),]
ibm_demo_rms[which(ibm_demo_rms$site == "NewMexico"),"site"] <- "1New Mexico"
ibm_demo_rms[which(ibm_demo_rms$site == "Arizona"),"site"] <- "2Arizona"
ibm_demo_rms[which(ibm_demo_rms$site == "Kansas"),"site"] <- "3Kansas"
ibm_demo_rms[which(ibm_demo_rms$site == "Montana"),"site"] <- "4Montana"
ibm_demo_rms[which(ibm_demo_rms$site == "Idaho"),"site"] <- "5Idaho"

ggplot(ibm_demo_rms, aes(x=(expansion^2), y=avg_synch, color=site))+
  geom_hline(data=theoretical_preds_and_obs, aes(yintercept=envonly_pred), color="grey15", linetype=3)+
  geom_hline(data=theoretical_preds_and_obs, aes(yintercept=demonly_pred), color="grey15", linetype=2)+
  geom_line(aes(group=experiment))+
  geom_point(aes(shape=experiment), size=3)+
  geom_errorbar(aes(ymin=lo_synch, ymax=up_synch), width=1.5)+
  facet_wrap("site", nrow=1)+
  scale_color_manual(values=site_colors, labels=site_labels, name="")+
  xlab(expression(paste("Simulated Area (", m^2,")")))+
  ylab("Synchrony of Species' Growth Rates")+
  scale_y_continuous(limits=c(0,1))+
  scale_shape_discrete(name="",labels=c("D.S. Only", "D.S + E.S."))+
  # scale_linetype_discrete(name="",labels=c("No E.S", "All Drivers"))+
  guides(color=FALSE)+
  theme_few()+
  theme(legend.position=c(0.1,0.2),
        legend.background = element_rect(fill = NA))

ggsave("../docs/components/ibm_sims_across_landscape.png", width = 12, height = 3, units="in", dpi=75)




####
####  SUPPLEMENTAL FIGURES 2-3
####
##  Make supplement plots for percent cover ------------
ipm_synch_all <- ddply(synch_df, .(site, experiment, typesynch), summarise,
                       mean_synch = mean(synch),
                       up_synch = quantile(synch, 0.95),
                       lo_synch = quantile(synch, 0.05))
ipm_synch <- subset(ipm_synch_all, typesynch=="abund_synch")

##  Read in IBM results -----
ibm_synch_all <- readRDS("../results/ibm_sims_collated.RDS") 
ibm_synch_agg <- ddply(ibm_synch_all, .(experiment, site, expansion, typesynch), summarise,
                       avg_synch = mean(synch),
                       up_synch = quantile(synch, 0.95),
                       lo_synch = quantile(synch, 0.05))
ibm_synch <- subset(ibm_synch_agg, expansion==5 & typesynch=="Cover (%)")

##  Combine results -----
sim_names <- c("All Drivers", "No D.S.", "No E.S.", "No Comp.", "No Comp. + No D.S.", "No Comp. + No E.S.")
sim_names_order <- paste0(c(3,1,5,4,2,6),sim_names)
site_names <- unique(ipm_synch$site)
nsites <- length(site_names)
site_labels <- site_names
site_labels[which(site_labels=="NewMexico")] <- "New Mexico"
site_labels_order <- paste0(c(2,5,3,4,1), site_labels)

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

plot_df <- data.frame(site = rep(site_labels_order, times=length(sim_names)),
                      simulation = rep(sim_names_order, each=nsites),
                      synchrony = all_experiments,
                      upper_synch = all_ups,
                      lower_synch = all_downs)


site_labels <- site_labels[order(site_labels_order)]
sim_labels <- sim_names[order(sim_names_order)]

# Bring in theoretical predictions and observed values
theory_env <- readRDS("../results/envstoch_predictions.RDS")
theory_env <- theory_env[2:nrow(theory_env), c("site", "cover_prediction")]
metrics <- readRDS("../results/calculate_site_metrics.RDS")
theory_demo <- metrics[,c("demo_only_exp_coversynch", "obs_cover_synch", "site")]
theory_demo[which(theory_demo$site=="New Mexico"),"site"] <- "NewMexico"
theoretical_preds_and_obs <- merge(theory_env, theory_demo)
colnames(theoretical_preds_and_obs) <- c("site", "envonly_pred", "demonly_pred", "obs")
theoretical_preds_and_obs$site <- site_labels_order

##  Make the main (all sims) plot -----
pointocols <- rep("grey25",3)
ggplot(data=plot_df)+
  geom_bar(data=plot_df, aes(x=simulation, y=synchrony, fill=site, color=site),
           stat="identity")+
  geom_errorbar(data=plot_df, aes(x=simulation, y=synchrony, fill=site, color=site, 
                                  ymin=lower_synch, ymax=upper_synch), 
                color="white", size=1, width=0.25)+
  geom_errorbar(data=plot_df, aes(x=simulation, y=synchrony, fill=site, color=site, 
                                  ymin=lower_synch, ymax=upper_synch), width=0.25)+
  geom_point(data=theoretical_preds_and_obs, aes(x=3, y=obs), size=3.5, color="white", shape=15)+
  geom_point(data=theoretical_preds_and_obs, aes(x=3, y=obs), size=3, color=pointocols[1], shape=15)+
  geom_point(data=theoretical_preds_and_obs, aes(x=2, y=envonly_pred), size=3.5, color="white", shape=16)+
  geom_point(data=theoretical_preds_and_obs, aes(x=2, y=envonly_pred), size=3, color=pointocols[2], shape=16)+
  geom_point(data=theoretical_preds_and_obs, aes(x=6, y=demonly_pred), size=3.5, color="white", shape=17)+
  geom_point(data=theoretical_preds_and_obs, aes(x=6, y=demonly_pred), size=3, color=pointocols[3], shape=17)+
  facet_wrap("site", nrow=1)+
  scale_fill_manual(values=site_colors, labels=site_labels, name="")+
  scale_color_manual(values=site_colors, labels=site_labels, name="")+
  xlab("Simulation Experiment")+
  ylab("Synchrony of Species' Percent Cover")+
  scale_x_discrete(labels=sim_labels)+
  scale_y_continuous(limits=c(0,1))+
  theme_few()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=FALSE,color=FALSE)

ggsave("../docs/components/figureS2_simscover.png", width = 10, height = 4, units="in", dpi=75)


##  Make demographic stochasiticty plot for all landscape sizes -----
ibm_demo_rms <- subset(ibm_synch_agg, typesynch=="Cover (%)")
ibm_demo_rms <- ibm_demo_rms[which(ibm_demo_rms$experiment %in% c("fluctnointer", "constnointer")),]
ibm_demo_rms[which(ibm_demo_rms$site == "NewMexico"),"site"] <- "1New Mexico"
ibm_demo_rms[which(ibm_demo_rms$site == "Arizona"),"site"] <- "2Arizona"
ibm_demo_rms[which(ibm_demo_rms$site == "Kansas"),"site"] <- "3Kansas"
ibm_demo_rms[which(ibm_demo_rms$site == "Montana"),"site"] <- "4Montana"
ibm_demo_rms[which(ibm_demo_rms$site == "Idaho"),"site"] <- "5Idaho"

ggplot(ibm_demo_rms, aes(x=(expansion^2), y=avg_synch, color=site))+
  geom_hline(data=theoretical_preds_and_obs, aes(yintercept=envonly_pred), color="grey15", linetype=3)+
  geom_hline(data=theoretical_preds_and_obs, aes(yintercept=demonly_pred), color="grey15", linetype=2)+
  geom_line(aes(group=experiment))+
  geom_point(aes(shape=experiment), size=3)+
  geom_errorbar(aes(ymin=lo_synch, ymax=up_synch), width=2)+
  facet_wrap("site", nrow=1)+
  scale_color_manual(values=site_colors, labels=site_labels, name="")+
  xlab(expression(paste("Simulated Area (", m^2,")")))+
  ylab("Synchrony of Species' Percent Cover")+
  scale_y_continuous(limits=c(0,1))+
  scale_shape_discrete(name="",labels=c("D.S. Only", "D.S + E.S."))+
  # scale_linetype_discrete(name="",labels=c("No E.S", "All Drivers"))+
  guides(color=FALSE)+
  theme_few()+
  theme(legend.position=c(0.1,0.2),
        legend.background = element_rect(fill = NA))

ggsave("../docs/components/figureS3_ibmcover.png", width = 10, height = 3, units="in", dpi=75)




####
####  VARIABILITY v. SYNCHRONY PLOTS -- ALL SIMULATION EXPERIMENTS (Fig. S4)
####

##  Read in IPM results ---------
output_list <- readRDS("../results/ipm_comp_nocomp_sims.RDS")
mlist <- melt(output_list)
colnames(mlist)[1:3] <- c("year", "species", "cover")
sites <- unique(mlist$L1)
sims <- unique(mlist$L2)
synch_df <- list() # empty storage
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
      tmpsynch <- get_ipm_synchrony(tmp)
      tmp_pgr_synch <- as.numeric(tmpsynch$pgr_synchrony["obs"])
      tmpcast <- dcast(tmp, year~species, value.var = "cover")
      tmp_abund_synch <- as.numeric(community.sync(tmpcast[2:ncol(tmpcast)])[1])
      tmp_abund_cv <- sd(rowSums(tmpcast[2:ncol(tmpcast)])) / mean(rowSums(tmpcast[2:ncol(tmpcast)]))
      tmp_df <- data.frame(site=dosite, experiment=dosim, 
                           bootnum=i, pgr_synch=tmp_pgr_synch, 
                           abund_synch=tmp_abund_synch,
                           abund_cv=tmp_abund_cv)
      synch_df <- rbind(synch_df, tmp_df)
    }# end boots loop
  }# end experiment/sim loop
}# end site loop


##  Read in IBM results ---
totSims <- 100   # number of simulations per site 
totT <- 100      # time steps of simulation
burn.in <- 25    # time steps to discard before calculating cover values

### Function for inverse of %in%
"%w/o%" <- function(x, y) x[!x %in% y] # x without y

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


output_df <- list()

for(do_exp in unique(files_df$experiment)){
  exp_files <- files_df[which(files_df$experiment==do_exp),"filename"]
  
  for(do_site in site_names){ # loop over sites
    tmp_files <- exp_files[grep(do_site, exp_files)]
    nsims <- length(tmp_files)
    
    for(i in 5){ # just do largest plot size, 5
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
          cv.abund <- sd(rowSums(tmpsynch[burn.in:totT,colids])) / mean(rowSums(tmpsynch[burn.in:totT,colids]))
          synch.pgr <- as.numeric(community.sync(obs_gr)[1])
          tmpD <- data.frame(experiment=do_exp, site=do_site, expansion=i, 
                             run=count, pgr_synch=synch.pgr, 
                             abund_synch=synch.abund, abund_cv=cv.abund)
          output_df <- rbind(output_df, tmpD)
          count <- count+1
        } # end loop over coexistence runs
        
      } # end if/then for coexistence runs
      
    }# end plot size sim loop
    
  }# end site loop
  
}# end experiment loop


### Combine IPM and IBM results
plot_cv_df <- rbind(output_df[,c("site", "experiment", "pgr_synch", "abund_synch", "abund_cv")],
                    synch_df[,c("site", "experiment", "pgr_synch", "abund_synch", "abund_cv")])

site_names <- as.character(unique(plot_cv_df$site))
nsites <- length(site_names)
site_labels <- site_names
site_labels[which(site_labels=="NewMexico")] <- "New Mexico"
site_labels_order <- paste0(c(2,5,3,4,1), site_labels)
site_labels <- site_labels[order(site_labels_order)]

site_order_nums <- c(2,5,3,4,1)
plot_cv_df$sites_order <- NA
for(i in 1:length(site_names)){
  do_site <- sort(site_names)[i]
  plot_cv_df[which(plot_cv_df$site==do_site), "sites_order"] <- paste0(site_order_nums[i],do_site)
}

exp_labs <- c("No E.S.", "No Comp. + No E.S.", "No D.S.", "No Comp. + No D.S.", "All Drivers", "No Comp.")


exp_colors <- c("#c7533b",
                "#8e4cbe",
                "#66a659",
                "#ae517a",
                "#9c7e40",
                "#6984b3")

library(gridExtra)

g1 <- ggplot(plot_cv_df, aes(x=abund_synch, y=abund_cv, color=experiment))+
  geom_point(shape=19, alpha=0.4)+
  geom_point(shape=1, alpha=0.6)+
  stat_smooth(method="lm", se=FALSE, size=1, linetype=1)+
  facet_wrap("sites_order", nrow=1)+
  scale_x_continuous(limits=c(0,1), labels=c("0","0.25","0.50","0.75","1"))+
  scale_color_manual(values=exp_colors, name="", labels=exp_labs)+
  # scale_fill_manual(values=c("lightcoral", "skyblue"), name="")+
  xlab("Synchrony of percent cover")+
  ylab("CV of total community cover")+
  theme_few()+
  guides(fill=FALSE)+
  theme(legend.key=element_rect(size=1),
        legend.key.size = unit(0.8, "lines"))

g2 <- ggplot(plot_cv_df, aes(x=pgr_synch, y=abund_cv, color=experiment))+
  geom_point(shape=19, alpha=0.4)+
  geom_point(shape=1, alpha=0.6)+
  stat_smooth(method="lm", se=FALSE, size=1, linetype=1)+
  facet_wrap("sites_order", nrow=1)+
  scale_x_continuous(limits=c(0,1), labels=c("0","0.25","0.50","0.75","1"))+
  scale_color_manual(values=exp_colors, name="", labels=exp_labs)+
  # scale_fill_manual(values=c("lightcoral", "skyblue"), name="")+
  xlab("Synchrony of per capita growth rates")+
  ylab("CV of total community cover")+
  theme_few()+
  guides(fill=FALSE)+
  theme(legend.key=element_rect(size=1),
        legend.key.size = unit(0.8, "lines"))

setwd("~/Repos/community_synchrony/docs/components/")
png("cv_vs_synchrony.png", height=5, width=12, units = "in", res = 72)
grid.arrange(g1, g2, nrow=2)
dev.off()



### Plot average CV versus average synchrony
# cv_synch_agg <- ddply(plot_cv_df, .(sites_order,experiment), summarise,
#                       avg_pgr_synch = mean(pgr_synch),
#                       avg_abund_synch = mean(abund_synch),
#                       avg_abund_cv = mean(abund_cv))
# 
# ggplot(cv_synch_agg, aes(x=avg_abund_synch, y=avg_abund_cv, color=experiment))+
#   geom_point(shape=19, alpha=0.4)+
#   geom_point(shape=1, alpha=0.6)+
#   facet_wrap("sites_order", nrow=1)+
#   scale_x_continuous(limits=c(0,1), labels=c("0","0.25","0.50","0.75","1"))+
#   scale_color_manual(values=exp_colors, name="", labels=exp_labs)+
#   # scale_fill_manual(values=c("lightcoral", "skyblue"), name="")+
#   xlab("Synchrony of per capita growth rates")+
#   ylab("CV of total community cover")+
#   theme_few()+
#   guides(fill=FALSE)+
#   theme(legend.key=element_rect(size=1),
#         legend.key.size = unit(0.8, "lines"))
