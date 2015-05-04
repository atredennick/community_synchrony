##  Script to source recruitment regression function

rm(list=ls(all=TRUE))
library(communitySynchrony)

# Set up some paths and lists
path_to_data <- "../data/"
site_list <- list.files(path_to_data)
removes <- c(grep("*.csv", site_list),
             grep("Climate", site_list),
             grep("Crowding", site_list))
site_list <- site_list[-removes]

for(do_site in site_list){
  species_list <- list.files(paste(path_to_data, do_site, "/", sep=""))
  
  i <- 1 #counter
  for(do_species in species_list){
    tmpfile <- paste(path_to_data, do_site, "/", do_species,"/recArea.csv",sep="")
    tmpD <- read.csv(tmpfile)
    
    # Add group info for each site individually, as needed
    if(do_site=="Arizona")
      tmpD$Group=as.factor(substr(tmpD$quad,1,1))
    if(do_site=="Kansas")
      tmpD$Group=as.numeric(tmpD$Group)-1
    if(do_site=="Montana")
      tmpD$Group=as.factor(substr(tmpD$quad,1,1)) 
    if(do_site=="NewMexico")
      tmpD$Group=as.factor(substr(tmpD$quad,1,1))
   
    tmpD$Group <- as.factor(tmpD$Group)
    tmpD <- tmpD[,c("quad","year","NRquad","totParea","Group")]
    names(tmpD)[3] <- paste("R.",species_list[i],sep="")
    names(tmpD)[4] <- paste("cov.",species_list[i],sep="")
    if(i==1){
      D=tmpD
    }else{
      D=merge(D,tmpD,all=T)
    }
    i <- i+1
  } # end species loop
  D[is.na(D)] <- 0  # replace missing values
  
  
} # end site loop

