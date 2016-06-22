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
rec_param_summaries <- readRDS("../results/recruit_parameters_yrlycomp.RDS")
sites <- names(rec_param_summaries)
plot_list <- list()
for(do_site in sites){
  site_summary <- rec_param_summaries[[do_site]]
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