##  Script to source vital rate regression functions

library(communitySynchrony)

# Loop through sites and then species

path_to_data <- "../data/"
site_list <- list.files(path_to_data)
removes <- c(grep("*.csv", site_list),
             grep("Climate", site_list),
             grep("Crowding", site_list))
site_list <- site_list[-removes]



growth_params <- list()
for(spp in sppList){
  tmp <- get_growth_params(dataframe = all_dataframes[[spp]],
                           crowd_mat = all_crowding[[spp]])
  growth_params[[spp]] <- tmp
}
