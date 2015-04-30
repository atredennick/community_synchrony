##  Script to call crowding functions and save results
library(communitySynchrony)

# Change this if directory structure is different
path_to_data <- "../data/"

site_list <- list.files("../data/")
removes <- c(grep("*.csv", site_list),
             grep("Climate", site_list),
             grep("Crowding", site_list))
site_list <- site_list[-removes]

growth_alphas <- read.csv(paste(path_to_data, "alpha_list_growth.csv", sep=""))
survival_alphas <- read.csv(paste(path_to_data, "alpha_list_survival.csv", sep=""))

crowd_list <- list()

# Loop through sites and get crowding for growth and survival
#   save each output in the respective data folder
for(do_site in site_list){
  tmp_alpha_grow <- subset(growth_alphas, Site==do_site)$Alpha
  tmp_alpha_surv <- subset(survival_alphas, Site==do_site)$Alpha
  
  tmp_crowd_grow <- estimate_crowding(site = do_site, 
                                      data_path = path_to_data,
                                      alphas = tmp_alpha_grow,
                                      vital_rate = "growth")
  
  tmp_crowd_surv <- estimate_crowding(site = do_site, 
                                      data_path = path_to_data,
                                      alphas = tmp_alpha_surv,
                                      vital_rate = "survival")
  
  # Save output
  grow_list[[do_site]] <- tmp_crowd_grow
  surv_list[[do_site]] <- tmp_crowd_surv
}

names(grow_list) <- site_list
names(surv_list) <- site_list

saveRDS(grow_list, "../data/Crowding/crowd_grow_list.RDS")
saveRDS(surv_list, "../data/Crowding/crowd_survival_list.RDS")