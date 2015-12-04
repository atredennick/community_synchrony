##  Devel script for plotting/showing correlations
##  among spp random year effects. // NOW COMMENTED OUT

##  Devel script for sychrony of yearly intrinsic growth rates

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 10-29-2015


####
####  Clean workspace; Load libraries ------------------------------------------
####
rm(list = ls()) # wipe the workspace clean
library(reshape2)
library(plyr)
library(ggplot2)
library(synchrony)



####
####  Read in growth rate matrices ---------------------------------------------
####
pgr_list <- readRDS("../results/ipm_yearly_pgr.RDS")
site_names <- names(pgr_list)
out_df <- data.frame(site=NA, pgr_synch=NA)
for(do_site in site_names){
  tmp_pgrs <- pgr_list[[do_site]]
  tmp_synch <- as.numeric(community.sync(tmp_pgrs)[1])
  tmp_df <- data.frame(site=do_site, pgr_synch=tmp_synch)
  out_df <- rbind(out_df, tmp_df)
}
pgr_synch <- out_df[2:nrow(out_df),]



####
####  Read in regression fits --------------------------------------------------
####
all_files <- list.files("../results/")
param_files <- all_files[grep("param*", all_files)]
num_vitals <- length(param_files)

list_yr_effects <- list()
for(j in 1:num_vitals){
  param_list <- readRDS(paste0("../results/",param_files[j]))
  vital_name <- unlist(strsplit(param_files[j], "_"))[1]
  
  if(vital_name != "recruit"){
    
    site_names <- names(param_list)
    site_year_list <- list()
    for(do_site in site_names){
      tmp_site <- param_list[[do_site]]
      num_spp <- length(tmp_site) 
      tmp_yr_mat <- matrix(ncol = num_spp, nrow=nrow(tmp_site[[1]]))
      for(i in 1:num_spp){
        spp_df <- tmp_site[[i]]
        tmp_yr_mat[,i] <- spp_df$Intercept.yr
      } # end species within site loop
      site_year_list[[do_site]] <- tmp_yr_mat
    } # end site within vital rate loop
    
  } # end NOT recruit loops
 
  if(vital_name == "recruit"){
    
    site_names <- names(param_list)
    site_year_list <- list()
    for(do_site in site_names){
      tmp_site <- param_list[[do_site]]
      tmp_yrs <- tmp_site[grep("intcpt.yr", rownames(tmp_site)),"Mean"]
      stringlist <- strsplit(names(tmp_yrs), ",")
      num_spp <- length(unique(sapply(stringlist, "[[", 2)))
      mat_yrs <- t(matrix(tmp_yrs, nrow=num_spp, byrow = TRUE))
      site_year_list[[do_site]] <- mat_yrs
    } # end site within vital rate loop
    
  } # end YES recruit loops
  
  list_yr_effects[[vital_name]] <- site_year_list
} # end vital rate loop


####
####  Calculate average correlation of year effects ----------------------------
####
vital_rates <- names(list_yr_effects)
out_df <- data.frame(site=NA, vital_rate=NA, corr=NA)
for(do_vital in vital_rates){
  tmp_vital <- list_yr_effects[[do_vital]]
  site_names <- names(tmp_vital)
  for(do_site in site_names){
    tmp_site <- tmp_vital[[do_site]]
    tmp_cor <- cor(tmp_site)
    avg_cor <- mean(tmp_cor[upper.tri(tmp_cor)])
    tmp_out <- data.frame(site=do_site, vital_rate=do_vital, corr=avg_cor)
    out_df <- rbind(out_df, tmp_out)
  }
}
cor_df <- out_df[2:nrow(out_df),]
cor_cast <- dcast(cor_df, site~vital_rate)


####
####  Calculate synchrony of random year effects -------------------------------
####
vital_rates <- names(list_yr_effects)
out_df <- data.frame(site=NA, vital_rate=NA, synch=NA)
for(do_vital in vital_rates){
  tmp_vital <- list_yr_effects[[do_vital]]
  site_names <- names(tmp_vital)
  for(do_site in site_names){
    tmp_site <- tmp_vital[[do_site]]
    site_synch <- as.numeric(community.sync(tmp_site)[1])
    tmp_out <- data.frame(site=do_site, vital_rate=do_vital, synch=site_synch)
    out_df <- rbind(out_df, tmp_out)
  }
}

out_df <- out_df[2:nrow(out_df),]
out_cast <- dcast(out_df, site~vital_rate)
rm(list=setdiff(ls(), c("out_cast", "out_df", "cor_df", 
                        "cor_cast", "pgr_synch")))


####
####  Save output for tables
####
outlist <- list(synch_long=out_df, synch_wide=out_cast,
                corr_long=cor_df, corr_wide=cor_cast,
                pgr_synch=pgr_synch)
# outfile = "path/here/file.RDS"
# saveRDS(outlist, outfile)



####
####  Plot pgr synch versus vital rate synchs
####
pgr_vr <- data.frame(site=outlist$synch_wide$site,
                     pgr_synch=outlist$pgr_synch$pgr_synch,
                     Survival=outlist$synch_wide$surv,
                     Growth=outlist$synch_wide$growth,
                     Recruitment=outlist$synch_wide$recruit)
pgr_vr_long <- melt(pgr_vr, id.vars = c("site","pgr_synch"))

ggplot(pgr_vr_long, aes(x=pgr_synch, y=value, group=variable))+
  geom_point()+
  geom_text(aes(label=site, y=value+0.02), size=3)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  facet_wrap("variable")+
  scale_x_continuous(limits=c(0,1))+
  scale_y_continuous(limits=c(0,1))+
  xlab("Synchrony of Yearly Per Capita Growth Rates")+
  ylab("Synchrony of Random Year Effects")+
  theme_bw()




