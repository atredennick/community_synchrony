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
  if(length(grep("_", species_list)) > 0){
    rms <- grep("_", species_list)
    species_list <- species_list[-rms]
  }
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
    
    # Add group info for each site individually, as needed
    if(do_site=="Arizona")
      D$Group=as.factor(substr(D$quad,1,1))
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
    
    W <- as.data.frame(W)
    W$xID <- D$X
    all_crowding[[spp]] <- W
  } # end species loop
  names(all_crowding) <- species_list
  return(all_crowding)
}