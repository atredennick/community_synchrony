##  Script to source vital rate regression functions

rm(list=ls(all=TRUE))
library(communitySynchrony)

# Set up some paths and lists
path_to_data <- "../data/"
site_list <- list.files(path_to_data)
removes <- c(grep("*.csv", site_list),
             grep("Climate", site_list),
             grep("Crowding", site_list))
site_list <- site_list[-removes]

# Load all crowding data list
crowd_growth <- readRDS(paste(path_to_data, "Crowding/crowd_grow_list.RDS", sep=""))
crowd_surv <- readRDS(paste(path_to_data, "Crowding/crowd_survival_list.RDS", sep=""))

# Load all alphas
grow_alphas <- read.csv(paste(path_to_data, "alpha_list_growth.csv", sep=""))
surv_alphas <- read.csv(paste(path_to_data, "alpha_list_survival.csv", sep=""))

# Loop through sites and then species
growth_params_biglist <- list()
for(do_site in site_list){
  species_list <- list.files(paste(path_to_data, do_site, "/", sep=""))
  alpha_now <- subset(grow_alphas, Site==do_site)
  alpha_now <- as.numeric(alpha_now$Alpha)
  
  growth_params <- list()
  for(do_species in species_list){
    growDfile <- paste(path_to_data, do_site, "/", do_species,"/growDnoNA.csv",sep="")
    growD <- read.csv(growDfile)
    #TODO -- check with Peter about this allEdge subset
    D <- growD  #subset(growD,allEdge==0)
    D$logarea.t0 <- log(D$area.t0)
    D$logarea.t1 <- log(D$area.t1)
    D$quad <- as.character(D$quad)
    
    # Add group info for each site individually, as needed
    if(do_site=="Arizona")
      D$Group=as.factor(substr(D$quad,1,1))
    if(do_site=="Kansas")
      D$Group=as.numeric(D$Group)-1
    if(do_site=="Montana")
      D$Group=as.factor(substr(D$quad,1,1)) 
    if(do_site=="NewMexico")
      D$Group=as.factor(substr(D$quad,1,1))
    
    # Get the years right for Kansas
    if(do_site=="Kansas")
      D <- subset(D, year<68)
    
    # Get correct crowding matrix
    crowd_growth_now <- crowd_growth[[do_site]][[do_species]]
    
    # Run through the function
    tmp <- get_growth_params(dataframe = D,
                             crowd_mat = crowd_growth_now,
                             alpha = alpha_now)
    
    # Save in temporary list
    growth_params[[do_species]] <- tmp
  } #end species loop
  # Save to the big list
  growth_params_biglist[[do_site]] <- growth_params
} #end site loop

