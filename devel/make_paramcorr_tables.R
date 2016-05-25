##  Script to create tables for parameter 
##  correlations from vital rate regressions

rm(list=ls())

library("xtable")

in_dir <- "../results/param_covariance/"
tmp_files <- list.files(in_dir)
all_files <- tmp_files[grep(".RDS", tmp_files)]

if(grep("surv_grow_paramcorrs", list.files("../results/")) != 0){
  file.remove("../results/surv_grow_paramcorrs.txt")
}

for(i in 1:length(all_files)){
  tmp_corrs <- readRDS(paste0(in_dir,all_files[i]))
  vital <- strsplit(all_files[i], "_")[[1]][1]
  site <- strsplit(all_files[i], "_")[[1]][2]
  species <- strsplit(all_files[i], "_")[[1]][3]
  cap <- paste(site,vital,species)
  cat(print(xtable(tmp_corrs, digits=-2, caption = cap), caption.placement = "top"), file = "../results/surv_grow_paramcorrs.txt", append = TRUE, sep=" ")
}


