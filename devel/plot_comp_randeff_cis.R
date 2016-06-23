################################################################################
##  plot_comp_randeff_cis.R: plots the 95% CIs for the yearly random offsets  ##
##  on competition for survival, growth, and recruitment.                     ##
################################################################################

rm(list=ls(all.names = TRUE)) # clear the workspace
library(ggplot2)
library(gridExtra)


####
####  SURVIVAL
####
surv_param_summaries <- readRDS("../results/surv_summaries_yrlycomp_list.RDS")
sites <- names(surv_param_summaries)
plot_list <- list()
for(do_site in sites){
  site_summary <- surv_param_summaries[[do_site]]
  species <- names(site_summary)
  for(do_species in species){
    species_summary <- site_summary[[do_species]]
    comp_ids <- grep("W", names(species_summary))
    for(i in comp_ids){
      comp_summary <- species_summary[[i]]
      colnames(comp_summary)[c(4,6)] <- c("low", "high")
      comp_summary$year <- 1:nrow(comp_summary)
      new_plot <- ggplot(data=comp_summary, aes(x=year, y=mean))+
        geom_hline(aes(yintercept=0), color="red")+
        geom_errorbar(aes(ymin=low, ymax=high), width=0.25)+
        # geom_point(size=4)+
        ylab("yearly competition effect")+
        ggtitle(paste(do_site,do_species,names(species_summary)[i]))+
        theme(plot.title = element_text(size = 12),
              axis.title = element_text(size=8))
      plot_list <- c(plot_list, list(new_plot))
    } # end competition coefficient loop
  } # end species loop
} # end site loop

png(filename = "randyear_compCIs_survival.png", width = 1280, height = 960, units = "px")
do.call(grid.arrange, c(plot_list, list(ncol = 7)))
dev.off()



####
####  GROWTH
####
grow_param_summaries <- readRDS("../results/growth_summaries_yrlycomp_list.RDS")
sites <- names(grow_param_summaries)
plot_list <- list()
for(do_site in sites){
  site_summary <- grow_param_summaries[[do_site]]
  species <- names(site_summary)
  for(do_species in species){
    species_summary <- site_summary[[do_species]]
    comp_ids <- grep("W", names(species_summary))
    for(i in comp_ids){
      comp_summary <- species_summary[[i]]
      colnames(comp_summary)[c(4,6)] <- c("low", "high")
      comp_summary$year <- 1:nrow(comp_summary)
      new_plot <- ggplot(data=comp_summary, aes(x=year, y=mean))+
        geom_hline(aes(yintercept=0), color="red")+
        geom_errorbar(aes(ymin=low, ymax=high), width=0.25)+
        # geom_point(size=4)+
        ylab("yearly competition effect")+
        ggtitle(paste(do_site,do_species,names(species_summary)[i]))+
        theme(plot.title = element_text(size = 12),
              axis.title = element_text(size=8))
      plot_list <- c(plot_list, list(new_plot))
    } # end competition coefficient loop
  } # end species loop
} # end site loop

png(filename = "randyear_compCIs_growth.png", width = 1280, height = 960, units = "px")
do.call(grid.arrange, c(plot_list, list(ncol = 7)))
dev.off()



####
####  RECRUITMENT
####
# density-dependence (dd) is a 3-D array: i focal species, j competing species, k year
rec_param_summaries <- readRDS("../results/recruit_quants_yrlycomp.RDS")
sites <- names(rec_param_summaries)
plot_list <- list()
for(do_site in sites){
  site_summary <- rec_param_summaries[[do_site]]
  comp_summary <- site_summary[grep("dd\\[",rownames(site_summary)),]
  species_numid <- unique(substr(rownames(comp_summary), 4, 4))
  species_names <- names(surv_param_summaries[[do_site]])
  
  for(i in 1:length(species_numid)){
    focal_species_rows <- grep(paste0("dd\\[", species_numid[i]), rownames(comp_summary))
    tmp_densdep <- as.data.frame(comp_summary[focal_species_rows,])
    colnames(tmp_densdep)[c(1,5)] <- c("low", "high")
    tmp_densdep$year <- c(1:nrow(tmp_densdep))
    tmp_densdep$focal <- substr(rownames(tmp_densdep), 4, 4)
    tmp_densdep$competitor <- substr(rownames(tmp_densdep), 6, 6)
    tmp_densdep$intra <- "yes"
    tmp_densdep[which(tmp_densdep$focal != tmp_densdep$competitor),"intra"] <- "no"
    new_plot <- ggplot(data=tmp_densdep, aes(x=year, color=intra))+
      geom_hline(aes(yintercept=0), color="red")+
      geom_errorbar(aes(ymin=low, ymax=high), width=0.25)+
      # geom_point(size=4)+
      ylab("yearly competition effect")+
      ggtitle(paste(do_site,species_names[i]))+
      theme(plot.title = element_text(size = 12),
            axis.title = element_text(size=8))
    plot_list <- c(plot_list, list(new_plot))
  } # end focal species loop
} # end site loop

png(filename = "randyear_compCIs_rec.png", width = 1280, height = 960, units = "px")
do.call(grid.arrange, c(plot_list, list(ncol = 4)))
dev.off()