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

####
#### GROWTH -------------------------------------------------------------------
####
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
    if(do_site=="Arizona"){
      D$Group=as.factor(substr(D$quad,1,1))
      D=subset(D,year>=17)
    }
    
    if(do_site=="Kansas")
      D$Group=as.numeric(D$Group)-1
    
    if(do_site=="Montana"){
      ##then we moved some specific points:
      tmp2<-which(D$quad=="A12" & D$year==44)
      tmp3<-which(D$quad=="B1"  & D$year==44)
      tmp41<-which(D$quad=="E4" & D$year==33) 
      tmp42<-which(D$quad=="E4" & D$year==34) 
      tmp43<-which(D$quad=="E4" & D$year==43)
      tmp44<-which(D$quad=="E4" & D$year==44)
      tmpONE<-c(tmp2,tmp3,tmp41,tmp42,tmp43,tmp44)
      if(length(tmpONE)>0) D<-D[-tmpONE,]
      D$Group=as.factor(substr(D$quad,1,1)) 
    }
    
    if(do_site=="NewMexico")
      D$Group=as.factor(substr(D$quad,1,1))
    
    # Get the years right for Kansas
    if(do_site=="Kansas"){
      D <- subset(D, year<68)
      ##to remove some points:
      #for q25
      tmp1<-which(D$quad=="q25")
      #for q27
      tmp2<-which(D$quad=="q27")
      #for q28
      tmp3<-which(D$quad=="q28")
      #for q30
      tmp4<-which(D$quad=="q30")
      #for q31
      tmp5<-which(D$quad=="q31" & (D$year<35 | D$year>39))
      #for q32
      tmp6<-which(D$quad=="q32" & (D$year<35 | D$year>41))
      tmp<-c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
      D<-D[-tmp,]
    }
      
    # Get correct crowding matrix
    crowd_growth_now <- crowd_growth[[do_site]][[do_species]]
    tmpmerge <- merge(D, crowd_growth_now, by.x="X", by.y="xID")
    tokeep <- grep("V", colnames(tmpmerge))
    crowd_growth_now <- as.matrix(tmpmerge[,tokeep])
    
    # Run through the function
    tmp <- get_growth_params_yrlcomp(dataframe = D,
                             crowd_mat = crowd_growth_now,
                             alpha = alpha_now)
    
    # Save in temporary list
    growth_params[[do_species]] <- tmp
  } #end species loop
  # Save to the big list
  growth_params_biglist[[do_site]] <- growth_params
} #end site loop

# Save the big parameter list
saveRDS(growth_params_biglist, "../results/growth_params_yrlycomp_list.RDS")



####
#### SURVIVAL -----------------------------------------------------------------
####
# Loop through sites and then species
surv_params_biglist <- list()
for(do_site in site_list){
  species_list <- list.files(paste(path_to_data, do_site, "/", sep=""))
  alpha_now <- subset(surv_alphas, Site==do_site)
  alpha_now <- as.numeric(alpha_now$Alpha)
  
  surv_params <- list()
  for(do_species in species_list){
    survDfile <- paste(path_to_data, do_site, "/", do_species,"/survD.csv",sep="")
    survD <- read.csv(survDfile)
    #TODO -- check with Peter about this allEdge subset
    D <- survD  #subset(growD,allEdge==0)
    D$logarea <- log(D$area)
    D$quad <- as.character(D$quad)
    
    # Add group info for each site individually, as needed
    if(do_site=="Arizona"){
      D$Group=as.factor(substr(D$quad,1,1))
      D=subset(D,year>=17)
    }
    if(do_site=="Kansas")
      D$Group=as.numeric(D$Group)-1
    if(do_site=="Montana"){
      ##then we moved some specific points:
      tmp2<-which(D$quad=="A12" & D$year==44)
      tmp3<-which(D$quad=="B1"  & D$year==44)
      tmp41<-which(D$quad=="E4" & D$year==33) 
      tmp42<-which(D$quad=="E4" & D$year==34) 
      tmp43<-which(D$quad=="E4" & D$year==43)
      tmp44<-which(D$quad=="E4" & D$year==44)
      D$Group=as.factor(substr(D$quad,1,1)) 
    }
    if(do_site=="NewMexico")
      D$Group=as.factor(substr(D$quad,1,1))
    
    # Get the years right for Kansas
    if(do_site=="Kansas"){
      D <- subset(D, year<68)
      ##to remove some points:
      #for q25
      tmp1<-which(D$quad=="q25")
      #for q27
      tmp2<-which(D$quad=="q27")
      #for q28
      tmp3<-which(D$quad=="q28")
      #for q30
      tmp4<-which(D$quad=="q30")
      #for q31
      tmp5<-which(D$quad=="q31" & (D$year<35 | D$year>39))
      #for q32
      tmp6<-which(D$quad=="q32" & (D$year<35 | D$year>41))
      tmp<-c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
      D<-D[-tmp,]
    }
    
    # Get correct crowding matrix
    crowd_surv_now <- crowd_surv[[do_site]][[do_species]]
    tmpmerge <- merge(D, crowd_surv_now, by.x="X", by.y="xID")
    tokeep <- grep("V", colnames(tmpmerge))
    crowd_surv_now <- as.matrix(tmpmerge[,tokeep])
    
    # Run through the function
    tmp <- get_survival_params_yrlcomp(dataframe = D,
                               crowd_mat = crowd_surv_now,
                               alpha = alpha_now)
    
    # Save in temporary list
    surv_params[[do_species]] <- tmp
  } #end species loop
  # Save to the big list
  surv_params_biglist[[do_site]] <- surv_params
} #end site loop

# Save the big parameter list
saveRDS(surv_params_biglist, "../results/surv_params_yrlycomp_list.RDS")

