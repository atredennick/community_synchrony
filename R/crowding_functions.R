##  Functions to estimate crowding effects on focal species

#' Estimate growth crowding effects on focal species
#' 
#' @author Andrew Tredennick
#' @param site Name of the focal site (character).
#' @param data_path The directory path leading to the data folder, 
#'        ending with "/" (character).
#' @param alphas The alpha values for each species at the site, in alphabetical
#'               order by species.
#' @param vital_rate Vital rate of interest. "growth" or "survival" 
#' @return List of crowding matrices, one list entry per species.

estimate_crowding <- function(site, data_path, alphas, vital_rate){
  path_to_files <- paste(data_path, site, "/", sep="")
  species_list <- sort(list.files(path_to_files))
  num_species <- length(species_list)
  alpha.effect <- alphas
  all_crowding <- list()
  
  if(vital_rate=="growth")
    filename <- "/growDnoNA.csv"
  if(vital_rate=="survival")
    filename <- "/survD.csv"
  
  #From here on, originally written by Peter Adler (so email him with your problems :))
  print("This may take a bit, but I'm working...I promise.")
  for(spp in 1:length(species_list)){
    doSpp=species_list[spp]
    Dfile=paste(path_to_files,doSpp,filename,sep="")
    D=read.csv(Dfile)
    D$quad <- as.character(D$quad)
    if(site=="Kansas")
      D <- subset(D, year<68)
    
    # remove outliers (large plants that obviously do not turn into tiny plants) for ARTR only
#     if(doSpp=="ARTR"){
#       tmp=which(D$quad=="Q23" & D$year==45 & D$trackID==67)
#       tmp=c(tmp,which(D$quad=="Q12" & D$year==55 & D$trackID==25))
#       tmp=c(tmp,which(D$quad=="Q26" & D$year==45 & D$trackID==73))
#       D=D[-tmp,]
#     }else{
#       D=D
#     }
    
    # calculate crowding 
    for(i in 1:length(species_list)){
      distDfile=paste(path_to_files,species_list[i],"/",species_list[i],"_genet_xy.csv",sep="")
      if(i==1){
        distD=read.csv(distDfile)
        distD$nbSpp=species_list[i]  
      }else{
        tmp=read.csv(distDfile)
        tmp$nbSpp=species_list[i] 
        distD=rbind(distD,tmp)
      }
    }
    
    distD=distD[,c("quad","year","trackID","area","nbSpp","x","y")]
    W=matrix(NA,dim(D)[1],length(species_list))

    for(i in 1:dim(D)[1]){
      tmpD=subset(distD,year==D$year[i] & quad==D$quad[i])
      focal=which(tmpD$trackID==D$trackID[i] & tmpD$nbSpp==doSpp)
      xx=tmpD$x[focal] ; yy=tmpD$y[focal]
      tmpD$distance=sqrt((xx-tmpD$x)^2+(yy-tmpD$y)^2)
      tmpD=subset(tmpD,distance>0)
      if(dim(tmpD)[1]>0){
        for(k in 1:length(species_list)){
          sppI=which(tmpD$nbSpp==species_list[k])
          if(length(sppI)>0){
            W[i,k]=sum(exp(-1*alpha.effect[k]*tmpD$distance[sppI]^2)*tmpD$area[sppI])         
          }else{
            W[i,k]=0
          }
        }
      }else{
        W[i,]=0
      }   
    }
    
    all_crowding[[spp]] <- W
  } # end species loop
  names(all_crowding) <- species_list
  return(all_crowding)
}