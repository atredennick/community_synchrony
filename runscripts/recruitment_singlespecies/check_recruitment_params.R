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
  nspp <- length(tmp_dds)
  tmp_ints <- tmp_params[grep("intcpt.mu", rownames(tmp_params)), "Mean"]
  
  tmp_df <- data.frame(site=rep(do_site,nspp),intercepts=tmp_ints, intra_crowd=tmp_dds)
  singlespp_params <- rbind(singlespp_params,tmp_df)
}


###  Plot the parameters against one another
png("../../docs/components/multi_single_rec_param_comparison.png", width = 7, height=4, units = "in", res = 120)
par(mfrow=c(1,2), tcl=-0.2, mgp=c(2,0.5,0))

plot(singlespp_params$intercepts, multispp_params$intercepts, pch=19,
     xlab="Estimate from single species model", ylab="Estimate from multi-species model",
     las=1,col=as.numeric(singlespp_params$site), main="Mean Intercept")
points(singlespp_params$intercepts, multispp_params$intercepts)
legend(0,6.2,legend = unique(singlespp_params$site), pt.bg = c(1:5), pch=21, cex = 0.5)
abline(0,1)

plot(singlespp_params$intra_crowd, multispp_params$intra_crowd, pch=19,
     xlab="Estimate from single species model", ylab="Estimate from multi-species model",
     las=1, col=as.numeric(singlespp_params$site), main="Intraspecific Crowding")
points(singlespp_params$intra_crowd, multispp_params$intra_crowd)
abline(0,1)

dev.off()
