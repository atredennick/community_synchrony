##  Script to look at differences of model parameters for recruitment regressions
##  with and without interspecific interactions

rm(list=ls())

site_names <- c("Arizona", "Idaho", "Kansas", "Montana", "NewMexico")
nsites <- length(site_names)

### Read in multi-species parameters
all_rec_params_multispp <- readRDS("../../results/recruit_parameters.RDS")
multispp_params <- list()

for(do_site in site_names){
  tmp_params <- all_rec_params_multispp[[do_site]]
  tmp_dds <- tmp_params[grep("dd", rownames(tmp_params)), "Mean"]
  nspp <- sqrt(length(tmp_dds))
  intra_dds <- diag(t(matrix(tmp_dds,nspp,nspp)))
  
  tmp_ints <- tmp_params[grep("intcpt.mu", rownames(tmp_params)), "Mean"]

  tmp_df <- data.frame(site=rep(do_site,nspp),intercepts=tmp_ints, intra_crowd=intra_dds)
  multispp_params <- rbind(multispp_params,tmp_df)
}

### Read in singles-species parameters
all_rec_params_multispp <- readRDS("recruit_ss_parameters.RDS")
singlespp_params <- list()

for(do_site in site_names){
  tmp_params <- all_rec_params_multispp[[do_site]]
  tmp_dds <- tmp_params[grep("dd", rownames(tmp_params)), "Mean"]
  nspp <- sqrt(length(tmp_dds))
  intra_dds <- diag(t(matrix(tmp_dds,nspp,nspp)))
  
  tmp_ints <- tmp_params[grep("intcpt.mu", rownames(tmp_params)), "Mean"]
  
  tmp_df <- data.frame(site=rep(do_site,nspp),intercepts=tmp_ints, intra_crowd=intra_dds)
  singlespp_params <- rbind(singlespp_params,tmp_df)
}