##  Functions that run all synchrony metrics given a dataset. The input
##    are the genet-level time series for a specific site. If quad-years
##    need to be removed, this should be done before feeding the data
##    to this function.

##  Outputs: -- time series of percent cover by species
##           -- time series of per capita growth rates by speciues
##           -- table of synchrony and stability metrics for the community

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Last update:  4.28.2015

#' Calculate synchrony and stability from genet-level time series
#' 
#' @param ts_data Time series dataframe of genet sizes. 'Bad' quad-years should already be removed.
#' @return A list of data frames for synchrony and stability metrics for the community.

get_comm_synchrony <- function(ts_data){
  ##  Get average size of all plants
  avg_cover <- mean(ts_data$totCover)
  
  ## aggregate the quadrat-level data for average observed cover
  # divide totCover by 100 to convert from m2 to percent cover in 1m2 plot
  ts_agg <- ddply(ts_data, .(year, species), summarise,
                  tot_cover = mean(totCover/100),
                  num_obs = length(totCover))
  ts_agg <- subset(ts_agg, num_obs>1) # ignore years where only 1 quadrat is observed
  
  ## Calculate stability
  # population stability (mean/sd)
  species_list <- unique(ts_agg$species)
  num_spp <- length(species_list)
  stability <- numeric(num_spp+1)
  obs_vector <- numeric(num_spp)
  for(i in 1:num_spp){ # loop through species for population stability
    tmp <- subset(ts_agg, species==species_list[i])
    stability[i] <- mean(tmp$tot_cover)/sd(tmp$tot_cover)
    obs_vector[i] <- nrow(tmp)
  } # end species looping for population stability
  
  # community stability (mean/sd of summed cover)
  # calculate total cover from the species-level data for each quadrat
  ts_sum <- ddply(ts_agg, .(year), summarise,
                  all_cover = sum(tot_cover))
  stability[num_spp+1] <- mean(ts_sum$all_cover)/sd(ts_sum$all_cover)
  stability <- as.data.frame(stability)
  stability$level <- c(species_list, "community")
  stability$num_observations <- c(obs_vector, nrow(ts_sum))
  
  
  ##  Calculate synchrony using Loreau and de Mazancourt metric
  # requires the R package 'synchrony'
  # Synchrony of abundance
  ts_mat <- dcast(ts_agg[,c("year", "species", "tot_cover")], 
                  formula = year~species, value.var="tot_cover")
  ts_mat <- ts_mat[which(complete.cases(ts_mat)==TRUE),]
  synch_abundance <- community.sync(ts_mat[2:(num_spp+1)]) #Loreau and de Mazancourt 2008 (Am Nat)
  
  # caclulate observed growth rates
  # create lagged data frame to only get observed yearly transitions
  lag_df <- ts_mat
  lag_df$lagyear <- lag_df$year+1
  colnames(lag_df)[2:(num_spp+1)] <- paste(colnames(lag_df)[2:(num_spp+1)],"_t0", sep="") 
  # merge the lag df with observed
  rm_col <- which(colnames(lag_df)=="year")
  merged_df <- merge(ts_mat, lag_df[,-rm_col], by.x = "year", by.y="lagyear")
  transitions <- nrow(merged_df)
  obs_gr <- matrix(nrow=transitions, ncol=num_spp)
  name_ids1 <- which(colnames(merged_df) %in% species_list)
  name_ids2 <- which(colnames(merged_df) %in% paste(species_list, "_t0", sep=""))
  for(i in 1:transitions){
    obs_gr[i,] <- as.numeric(log(merged_df[i,name_ids1]/merged_df[i,name_ids2]))
  }
  growth_rate_synchrony <- community.sync(obs_gr)
  
  obs_gr <- as.data.frame(obs_gr)
  colnames(obs_gr) <- species_list
  obs_gr$year <- merged_df$year
  
  ##  Total cover variability
  ts_tot <- apply(ts_mat[2:ncol(ts_mat)], 1, sum)
  comm_cv <- sd(ts_tot)/mean(ts_tot)
  
  
  ##  Expected synchrony under independent fluctuations
  # New formula from Claire
  ts_freq <- ddply(ts_data, .(year, species), summarise,
                   avg_totcov = mean(totCover))
  ts_freq_tmp <- dcast(ts_freq, year~species, value.var = "avg_totcov")
  ts_freq_tmp$total <- rowSums(ts_freq_tmp[,c(2:ncol(ts_freq_tmp))])
  ts_freq_wide <- colMeans(ts_freq_tmp[,c(2:(ncol(ts_freq_tmp)-1))] / ts_freq_tmp[,ncol(ts_freq_tmp)])
  expected_cover_synchrony <- 1 / (sum(ts_freq_wide^(1/2)))^2
  expected_pgr_synchrony <- sum(ts_freq_wide^-1) / (sum(ts_freq_wide^(-1/2)))^2
  
#   sigma <- numeric(num_spp)
#   sigma_sqr <- numeric(num_spp)
#   for(i in 1:num_spp){
#     sigma[i] <- sd(obs_gr[,i])
#     sigma_sqr[i] <- (sd(obs_gr[,i]))^2
#   }
#   expected_synchrony_ind_flucts <- (sum(sigma_sqr)) / ((sum(sigma))^2)
  
  
  ##  Output
  obs_gr <- melt(obs_gr, id.vars = "year")
  obs_gr <- obs_gr[with(obs_gr, order(year)), ]
  colnames(obs_gr) <- c("year","species","pgr")
  return(list(avg_size = avg_cover,
              stability = stability,
              variability = comm_cv,
              pgr_synchrony = growth_rate_synchrony,
              pgr_expected_synch_ind_flucts = expected_pgr_synchrony,
              cover_expected_synch_ind_flucts = expected_cover_synchrony,
              abund_synchrony = synch_abundance,
              growth_rates = obs_gr,
              percent_cover = ts_agg))
}




#' Calculate synchrony and stability from IPM time series
#' 
#' @param ts_data A melted dataframe of IPM cover values with three columns: year, species, cover.
#' @return A list of data frames for synchrony and stability metrics for the community.
get_ipm_synchrony <- function(ts_data){
  species_list <- unique(ts_data$species)
  num_spp <- length(species_list)
  
  ## Calculate stability
  # population stability (mean/sd)
  stability <- numeric(num_spp+1)
  obs_vector <- numeric(num_spp)
  for(i in 1:num_spp){ # loop through species for population stability
    tmp <- subset(ts_data, species==species_list[i])
    stability[i] <- mean(tmp$cover)/sd(tmp$cover)
    obs_vector[i] <- nrow(tmp)
  } # end species looping for population stability

  # community stability (mean/sd of summed cover)
  # calculate total cover from the species-level data for each quadrat
  ts_sum <- ddply(ts_data, .(year), summarise,
                  all_cover = sum(cover))
  stability[num_spp+1] <- mean(ts_sum$all_cover)/sd(ts_sum$all_cover)
  stability <- as.data.frame(stability)
  stability$level <- c(species_list, "community")
  stability$num_observations <- c(obs_vector, nrow(ts_sum))

  ##  Calculate synchrony using Loreau and de Mazancourt metric
  # requires the R package 'synchrony'
  # Synchrony of abundance
  ts_mat <- dcast(ts_data, formula = year~species, value.var="cover")
  ts_mat <- ts_mat[which(complete.cases(ts_mat)==TRUE),]
  synch_abundance <- community.sync(ts_mat[2:(num_spp+1)]) #Loreau and de Mazancourt 2008 (Am Nat)
  
  # caclulate observed growth rates
  # create lagged data frame to only get observed yearly transitions
  lag_df <- ts_mat
  lag_df$lagyear <- lag_df$year+1
  colnames(lag_df)[2:(num_spp+1)] <- paste(colnames(lag_df)[2:(num_spp+1)],"_t0", sep="") 
  # merge the lag df with observed
  rm_col <- which(colnames(lag_df)=="year")
  merged_df <- merge(ts_mat, lag_df[,-rm_col], by.x = "year", by.y="lagyear")
  transitions <- nrow(merged_df)
  obs_gr <- matrix(nrow=transitions, ncol=num_spp)
  name_ids1 <- which(colnames(merged_df) %in% species_list)
  name_ids2 <- which(colnames(merged_df) %in% paste(species_list, "_t0", sep=""))
  for(i in 1:transitions){
    obs_gr[i,] <- as.numeric(log(merged_df[i,name_ids1]/merged_df[i,name_ids2]))
  }
  growth_rate_synchrony <- community.sync(obs_gr)
  
  obs_gr <- as.data.frame(obs_gr)
  colnames(obs_gr) <- species_list
  obs_gr$year <- merged_df$year
  
  
  ##  Expected synchrony under independent fluctuations
  sigma <- numeric(num_spp)
  sigma_sqr <- numeric(num_spp)
  for(i in 1:num_spp){
    sigma[i] <- sd(obs_gr[,i])
    sigma_sqr[i] <- (sd(obs_gr[,i]))^2
  }
  expected_synchrony_ind_flucts <- (sum(sigma_sqr)) / ((sum(sigma))^2)
  
  
  ##  Output
  obs_gr <- melt(obs_gr, id.vars = "year")
  obs_gr <- obs_gr[with(obs_gr, order(year)), ]
  colnames(obs_gr) <- c("year","species","pgr")
  return(list(stability = stability,
              pgr_synchrony = growth_rate_synchrony,
              pgr_expected_synch_ind_flucts = expected_synchrony_ind_flucts,
              abund_synchrony = synch_abundance,
              growth_rates = obs_gr))
}
