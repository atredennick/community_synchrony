##  R script to calculate expected synchrony when only environmental
##  stochasticity is operating

## clear the workspace
rm(list=ls())
library(plyr)
library(reshape2)

##  Define functions
sync_biomass <- function(obs_spp_frequencies, mono_pairwise_spp_corrs) {
  pairwise_frequency <- expand.grid(obs_spp_frequencies, obs_spp_frequencies)
  corrs_vec <- as.numeric(mono_pairwise_spp_corrs)
  pairwise_frequency$Var3 <- corrs_vec
  sync <- sum(apply(pairwise_frequency, 1, prod))
  return(sync)
}

sync_growthrates <- function(mono_pairwise_spp_corrs) {
  sync <- sum(mono_pairwise_spp_corrs) / ncol(mono_pairwise_spp_corrs)^2
  return(sync)
}



####
####  Read in data and monoculture IPM results
####
ipm_results <- readRDS("../results/ipm_comp_nocomp_sims.RDS")
site_names <- c("Arizona", "Idaho", "Kansas", "Montana", "NewMexico")
synchrony_predictions <- data.frame(site=NA, 
                                    cover_prediction=NA, 
                                    growthrate_prediction=NA)
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
  if(do_site=="Kansas"){
    ks_data <- site_data
    tmp1<-which(ks_data$quad_data=="q25" & (ks_data$year<35 | ks_data$year>62))
    tmp2<-which(ks_data$quad_data=="q27")
    tmp3<-which(ks_data$quad=="q28")
    tmp4<-which(ks_data$quad=="q30")
    tmp5<-which(ks_data$quad=="q31" & (ks_data$year<35 | ks_data$year>39))
    tmp6<-which(ks_data$quad=="q32" & (ks_data$year<35 | ks_data$year>41))
    tmp<-c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
    site_data<-ks_data[-tmp,]
  }
  
  ts_freq <- ddply(site_data, .(year, species), summarise,
                   avg_totcov = mean(totCover))
  ts_freq_tmp <- dcast(ts_freq, year~species, value.var = "avg_totcov")
  ts_freq_tmp$total <- rowSums(ts_freq_tmp[,c(2:ncol(ts_freq_tmp))])
  obs_spp_frequencies <- colMeans(ts_freq_tmp[,c(2:(ncol(ts_freq_tmp)-1))] / ts_freq_tmp[,ncol(ts_freq_tmp)])
  
  site_ipms <- as.data.frame(ipm_results[[do_site]][["ENVNOINTER"]])
  site_ipms$year <- c(1:nrow(site_ipms))
  ts_mat <- site_ipms
  # caclulate observed growth rates
  # create lagged data frame to only get observed yearly transitions
  lag_df <- ts_mat
  lag_df$lagyear <- lag_df$year+1
  colnames(lag_df)[1:(num_spp)] <- paste(colnames(lag_df)[1:(num_spp)],"_t0", sep="") 
  # merge the lag df with observed
  rm_col <- which(colnames(lag_df)=="year")
  merged_df <- merge(ts_mat, lag_df[,-rm_col], by.x = "year", by.y="lagyear")
  transitions <- nrow(merged_df)
  obs_gr <- matrix(nrow=transitions, ncol=num_spp)
  name_ids1 <- which(colnames(merged_df) %in% spp_list)
  name_ids2 <- which(colnames(merged_df) %in% paste(spp_list, "_t0", sep=""))
  for(i in 1:transitions){
    obs_gr[i,] <- as.numeric(log(merged_df[i,name_ids1]/merged_df[i,name_ids2]))
  }
  mono_pairwise_spp_corrs <- as.matrix(cor(obs_gr))
  
  predicted_percentcover_synchrony <- sync_biomass(obs_spp_frequencies, mono_pairwise_spp_corrs)
  predicted_growthrate_synchrony <- sync_growthrates(mono_pairwise_spp_corrs)
  tmpdf <- data.frame(site=do_site,
                      cover_prediction=predicted_percentcover_synchrony,
                      growthrate_prediction=predicted_growthrate_synchrony)
  synchrony_predictions <- rbind(synchrony_predictions, tmpdf)
}

synchrony_predictions[2:nrow(synchrony_predictions),]
saveRDS(synchrony_predictions,"../results/envstoch_predictions.RDS")


