###########################################################################
##  calc_comp_year_corrs.R: script to calculate the correlation between  ##
##  year random effects on intercept, slope, and competitive effects.    ##
###########################################################################

rm(list=ls(all.names = TRUE)) # clear the workspace

library(reshape2)

####
####  GROWTH REGRESSIONS
####
all_growth_params <- readRDS("../results/growth_params_yrlycomp_list.RDS")
sites <- names(all_growth_params)

randeff_allsite_growth_cors <- list()
# do_site <- "Arizona"
for(do_site in sites){
  tmpsite_growth_params <- all_growth_params[[do_site]]
  species <- names(tmpsite_growth_params)
  
  randeff_out <- list()
  # do_species <- "BOER"
  for(do_species in species){
    tmpspecies_growth_params <- tmpsite_growth_params[[do_species]]
    year_columns <- grep(".yr", colnames(tmpspecies_growth_params))
    randeff_cors <- cor(tmpspecies_growth_params[,year_columns])
    row_torm <- grep("W1", rownames(randeff_cors))
    randeff_out[[do_species]] <- randeff_cors[-c(1,2,row_torm),]
  }
  randeff_allsite_growth_cors[[do_site]] <- randeff_out
}



####
####  SURVIVAL REGRESSIONS
####
all_surv_params <- readRDS("../results/surv_params_yrlycomp_list.RDS")
sites <- names(all_surv_params)

randeff_allsite_surv_cors <- list()
# do_site <- "Arizona"
for(do_site in sites){
  tmpsite_surv_params <- all_surv_params[[do_site]]
  species <- names(tmpsite_surv_params)
  
  randeff_out <- list()
  # do_species <- "BOER"
  for(do_species in species){
    tmpspecies_surv_params <- tmpsite_surv_params[[do_species]]
    year_columns <- grep(".yr", colnames(tmpspecies_surv_params))
    randeff_cors <- cor(tmpspecies_surv_params[,year_columns])
    row_torm <- grep("W1", rownames(randeff_cors))
    randeff_out[[do_species]] <- randeff_cors[-c(1,2,row_torm),]
  }
  randeff_allsite_surv_cors[[do_site]] <- randeff_out
}



####
####  RECRUITMENT REGRESSIONS
####
# density-dependence (dd) is a 3-D array: i focal species, j competing species, k year
# intercepts are 2-D: i year, j focal species
all_rec_params <- readRDS("../results/recruit_parameters_yrlycomp.RDS")
sites <- names(all_rec_params)

pdf("recruitment_intercept_comp_cors.pdf", onefile = TRUE, width=8, height=4)
randeff_allsite_rec_cors <- list()
for(do_site in sites){
  tmpsite_rec_params <- all_rec_params[[do_site]]
  densdep_rows <- grep("dd\\[", rownames(tmpsite_rec_params))
  densdep <- tmpsite_rec_params[densdep_rows, ]
  species_numid <- unique(substr(rownames(densdep), 4, 4))
  
  intercept_rows <- grep("intcpt.yr", rownames(tmpsite_rec_params))
  intercepts <- tmpsite_rec_params[intercept_rows, ]
  intercept_matrix <- as.data.frame(matrix(intercepts[,"Mean"], 
                                           nrow=nrow(intercepts)/length(species_numid),
                                           ncol=length(species_numid)))
  
  species_names <- names(all_surv_params[[do_site]])
  randeff_out <- list()
  for(i in 1:length(species_numid)){
    focal_species_rows <- grep(paste0("dd\\[", species_numid[i]), rownames(densdep))
    tmp_densdep <- densdep[focal_species_rows,]
    num_species <- length(species_numid)
    num_years <- nrow(tmp_densdep)/num_species
    tmp_densdep_matrix <- matrix(data = tmp_densdep[,"Mean"], 
                                 nrow=num_years, ncol=num_species, byrow = TRUE)
    tmp_densdep_matrix <- as.data.frame(tmp_densdep_matrix)
    colnames(tmp_densdep_matrix) <- paste0("W",species_numid,".yr")
    tmp_param_matrix <- cbind(intercept_matrix[,i], tmp_densdep_matrix)
    colnames(tmp_param_matrix)[1] <- "Intercept.yr"
    tmp_species_cors <- cor(tmp_param_matrix)
    melted <- melt(tmp_param_matrix, id.vars = "Intercept.yr")
    row_torm <- grep(paste0("W",i), rownames(tmp_species_cors))
    randeff_out[[species_names[i]]] <- tmp_species_cors[-c(1,row_torm),]
    
    gprint <- ggplot(melted, aes(x=Intercept.yr, y=value))+
      geom_point()+
      facet_wrap("variable", scales="free")+
      ggtitle(paste("recruitment", do_site, species_names[i]))
    print(gprint)
  }
  randeff_allsite_rec_cors[[do_site]] <- randeff_out
}
dev.off()


####
####  PLOT HISTOGRAMS OF CORRELATIONS BY REGRESSION
####
grow_max_cor <- 0
for(do_site in sites){
  tmp_grow <- randeff_allsite_growth_cors[[do_site]]
  species_names <- names(tmp_grow)
  for(do_species in species_names){
    tmpspp_grow <- tmp_grow[[do_species]]
    if(is.matrix(tmpspp_grow)==TRUE) { 
      if(max(tmpspp_grow[,1:2]) > grow_max_cor) grow_max_cor <- max(abs(tmpspp_grow[,1:2]))
    }
    if(is.matrix(tmpspp_grow)==FALSE) { 
      if(max(tmpspp_grow[1:2]) > grow_max_cor) grow_max_cor <- max(abs(tmpspp_grow[1:2]))
    }
  }
}

surv_max_cor <- 0
for(do_site in sites){
  tmp_surv <- randeff_allsite_surv_cors[[do_site]]
  species_names <- names(tmp_surv)
  for(do_species in species_names){
    tmpspp_surv <- tmp_surv[[do_species]]
    if(is.matrix(tmpspp_surv)==TRUE) { 
      if(max(tmpspp_surv[,1:2]) > surv_max_cor) surv_max_cor <- max(abs(tmpspp_surv[,1:2]))
    }
    if(is.matrix(tmpspp_surv)==FALSE) { 
      if(max(tmpspp_surv[1:2]) > surv_max_cor) surv_max_cor <- max(abs(tmpspp_surv[1:2]))
    }
  }
}

rec_max_cor <- 0
for(do_site in sites){
  tmp_rec <- randeff_allsite_rec_cors[[do_site]]
  species_names <- names(tmp_rec)
  for(do_species in species_names){
    tmpspp_rec <- tmp_rec[[do_species]]
    if(is.matrix(tmpspp_rec)==TRUE) { 
      if(max(tmpspp_rec[,1]) > rec_max_cor) rec_max_cor <- max(abs(tmpspp_rec[,1]))
    }
    if(is.matrix(tmpspp_rec)==FALSE) { 
      if(max(tmpspp_rec[1]) > rec_max_cor) rec_max_cor <- max(abs(tmpspp_rec[1]))
    }
  }
}


##  PRINT RESULTS
print(paste("Maximum absolute correlation for survival:", round(surv_max_cor,2)))
print(paste("Maximum absolute correlation for growth:", round(grow_max_cor,2)))
print(paste("Maximum absolute correlation for recruitment:", round(rec_max_cor,2)))



####
####  PLOT COMPETITION YEAR EFFECTS 
####
