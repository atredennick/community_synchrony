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
library(ggplot2)



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
# Growth and survival first
vital_mat_list <- list()
for(do_vital in vital_rates[c(1,2)]){
  vr_do <- vr_files[grep(do_vital,vr_files)]
  tmp_vital <- readRDS(paste0("../results/",vr_do))
  sites <- names(tmp_vital)
  mat_list <- list()
  for(do_site in sites){
    tmp <- tmp_vital[[do_site]]
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
  vital_mat_list[[do_vital]] <- mat_list 
}

# Recruitment
vr_do <- vr_files[grep("recruit",vr_files)]
tmp_vital <- readRDS(paste0("../results/",vr_do))
names_vital <- readRDS(paste0("../results/",vr_files[1]))
sites <- names(tmp_vital)
mat_list <- list()
for(do_site in sites){
  tmp <- tmp_vital[[do_site]]
  spp <- names(names_vital[[do_site]])
  tmpdd <- tmp[grep("dd", rownames(tmp)),"Mean"]
  tmpmat <- matrix(tmpdd,length(spp),length(spp), byrow = FALSE)
  rownames(tmpmat) <- spp
  colnames(tmpmat) <- spp
  mat_list[[do_site]] <- tmpmat
}
vital_mat_list[["recruit"]] <- mat_list 



####
####  Calculate Average Inter and Intraspecific Interactions -------------------
####
vital_comp_list <- list()
for(do_vital in vital_rates){
  tmp <- vital_mat_list[[do_vital]]
  sites <- names(tmp)
  comp_list <- list()
  for(do_site in sites){
    tmpsite <- tmp[[do_site]]
    avgintra <- mean(diag(tmpsite))
    avginter <- mean(c(tmpsite[upper.tri(tmpsite)], tmpsite[lower.tri(tmpsite)]))
    avgcomp <- matrix(c(avginter, avgintra),1,2)
    colnames(avgcomp) <- c("inter", "intra")
    comp_list[[do_site]] <- avgcomp
  }
  vital_comp_list[[do_vital]] <- comp_list
}



####
####  Plot Histogram of Spp Interactions ---------------------------------------
####
df <- data.frame(comp = as.numeric(unlist(vital_comp_list)),
                 comptype = rep(c("inter", "intra"),times=length(unlist(vital_comp_list))/2))
ggplot(df, aes(x=comp))+
  geom_bar(fill="grey25", color="white", binwidth=0.4)+
  facet_wrap("comptype", ncol=2)+
  theme_bw()



####
####  Plot Each Vital Rate Interaction Coefficients
####
df <- melt(unlist(vital_comp_list))
pieces <- strsplit(rownames(df), split = "[.]")
df$vitalrate <- sapply(pieces, "[", 1)
sites <- sapply(pieces, "[", 2)
df$site <- substr(sites, start = 1, stop = nchar(sites)-1)
df$interintra <- substr(sites, start = nchar(sites), stop = nchar(sites))
ggplot(df, aes(x=site, y=value, fill=interintra))+
  geom_bar(stat="identity", position="dodge")+
  facet_wrap("vitalrate", scales="free_y")



