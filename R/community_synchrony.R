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
  name_ids1 <- which(colnames(merged_df) %in% spp_list)
  name_ids2 <- which(colnames(merged_df) %in% paste(spp_list, "_t0", sep=""))
  for(i in 1:transitions){
    obs_gr[i,] <- as.numeric(log(merged_df[i,name_ids1]/merged_df[i,name_ids2]))
  }
  return(obs_gr)
}