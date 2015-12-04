##  Devel script making a table for the estimated
##  species interaction matrices at each site

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 12-04-2015


####
####  Clean Workspace an Load Libraries ----------------------------------------
####
rm(list = ls()) # wipe the workspace clean
library(reshape2)
library(plyr)



####
####  Get Vital Rate Satistical Results Files ----------------------------------
####
vital_rates <- c("growth", "surv", "recruit")
result_files <- list.files("../results/")
vr_files <-result_files[grep(paste(vital_rates,collapse="|"), 
                              list.files("../results/"))]



####
####  Loop Over Files and Extract Crowding Effects -----------------------------
####
test <- readRDS(paste0("../results/",vr_files[1]))
sites <- names(test)
mat_list <- list()
for(do_site in sites){
  tmp <- test[[do_site]]
  spp <- names(tmp)
  tmpmat <- matrix(0,nrow = length(spp), ncol = length(spp))
  for(i in 1:length(spp)){
    do_spp <- spp[i]
    tmpspp <- tmp[[do_spp]]
    tmpw <- tmpspp[1,grep("crowd*", colnames(tmpspp))]
    tmpmat[i,] <- as.numeric(tmpw)
  }
  rownames(tmpmat) <- spp
  colnames(tmpmat) <- spp
  mat_list[[do_site]] <- tmpmat
}

