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
library(ggplot2)
library(reshape2)
library(plyr)
library(synchrony)
library(communitySynchrony)


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
      tmp_yr_int <- matrix(ncol = num_spp, nrow=nrow(tmp_site[[1]]))
      tmp_yr_slope <- matrix(ncol = num_spp, nrow=nrow(tmp_site[[1]]))
      for(i in 1:num_spp){
        spp_df <- tmp_site[[i]]
        tmp_yr_int[,i] <- spp_df$Intercept.yr
        col2use <- grep("logarea", colnames(spp_df))
        col2use2 <- grep(".yr", colnames(spp_df[col2use]))
        tmp_yr_slope[,i] <- spp_df[,col2use[col2use2]]
      } # end species within site loop
      site_year_list[[do_site]][["intercept"]] <- tmp_yr_int
      site_year_list[[do_site]][["slope"]] <- tmp_yr_slope
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


# ## Plot histograms of year effects by vital rate
# vital_rates <- names(list_yr_effects)
# par(mfrow=c(3,5))
# for(do_vital in vital_rates){
#   tmp_vital <- list_yr_effects[[do_vital]]
#   site_names <- names(tmp_vital)
#   for(do_site in site_names){
#     tmp_site <- tmp_vital[[do_site]]
#     hist(tmp_site, main=paste(do_site,do_vital), xlab="Year Effect", xlim=c(-3,8))
#   }
# }

##  Get standard deviation of year effect posterior



####
####  Calculate average correlation of year effects ----------------------------
####
vital_rates <- names(list_yr_effects)
out_df <- data.frame(site=NA, vital_rate=NA, term=NA, corr=NA)
for(do_vital in vital_rates){
  tmp_vital <- list_yr_effects[[do_vital]]
  site_names <- names(tmp_vital)
  for(do_site in site_names){
    tmp_site <- tmp_vital[[do_site]]
    
    if(do_vital != "recruit"){
      for(i in 1:length(tmp_site)){
        tmp_cor <- cor(tmp_site[[i]])
        avg_cor <- mean(tmp_cor[upper.tri(tmp_cor)])
        tmp_out <- data.frame(site=do_site, vital_rate=do_vital, term=names(tmp_site)[i], corr=avg_cor)
        out_df <- rbind(out_df, tmp_out)
      }
    }
    
    if(do_vital == "recruit"){
      tmp_cor <- cor(tmp_site)
      avg_cor <- mean(tmp_cor[upper.tri(tmp_cor)])
      tmp_out <- data.frame(site=do_site, vital_rate=do_vital, term="intercept", corr=avg_cor)
      out_df <- rbind(out_df, tmp_out)
    }
    
  }
}
cor_df <- out_df[2:nrow(out_df),]
cor_cast <- dcast(cor_df, site+term~vital_rate)
# avg_env_resp <- apply(cor_cast[,3:5], MARGIN = 1, FUN = "mean")

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



####
####  Plot Synchrony of Yearly Per Capita Growth Rates vs. PGR Synchrony
####
output_list <- readRDS("../results/ipm_comp_nocomp_sims.RDS")
mlist <- melt(output_list)
colnames(mlist)[1:3] <- c("year", "species", "cover")
sites <- unique(mlist$L1)
sims <- unique(mlist$L2)
synch_df <- data.frame(site=NA, experiment=NA, bootnum=NA, 
                       pgr_synch=NA, abund_synch=NA)
boots <- 100
num_iters <- 13
for(dosite in sites){
  tmp_data <- subset(mlist, L1==dosite)
  for(dosim in sims){
    tmpsim <- subset(tmp_data, L2==dosim)[c("year", "species", "cover")]
    for(i in 1:boots){
      begin_year <- sample(x = 1:(max(tmpsim$year)-num_iters), 1)
      end_year <- begin_year+num_iters
      tmp <- subset(tmpsim, year %in% begin_year:end_year)
      tmpsynch <- get_ipm_synchrony(tmp)
      tmp_pgr_synch <- as.numeric(tmpsynch$pgr_synchrony["obs"])
      tmpcast <- dcast(tmp, year~species, value.var = "cover")
      tmp_abund_synch <- as.numeric(community.sync(tmpcast[2:ncol(tmpcast)])[1])
      tmp_df <- data.frame(site=dosite, experiment=dosim, 
                           bootnum=i, pgr_synch=tmp_pgr_synch, 
                           abund_synch=tmp_abund_synch)
      synch_df <- rbind(synch_df, tmp_df)
    }# end boots loop
  }# end experiment/sim loop
}# end site loop

synch_dftmp <- synch_df[2:nrow(synch_df),]
synch_df <- melt(synch_dftmp, id.vars = c("site", "experiment", "bootnum"))
colnames(synch_df) <- c("site", "experiment", "bootnum", "typesynch", "synch")

agg_synch <- ddply(synch_df, .(site, experiment, typesynch), summarise,
                   mean_synch = mean(synch),
                   up_synch = quantile(synch, 0.95),
                   lo_synch = quantile(synch, 0.05))

all_synch <- subset(agg_synch, experiment=="ENVINTER" & typesynch=="pgr_synch")

vr_synchs <- data.frame(Survival=outlist$synch_wide$surv,
                        Growth=outlist$synch_wide$growth,
                        Recruitment=outlist$synch_wide$recruit)
# vr_synchs <- data.frame(Survival=outlist$synch_wide$surv,
#                         Growth=outlist$synch_wide$growth)
vr_avg_synch <- rowMeans(vr_synchs)
all_synch$vravg <- vr_avg_synch
all_synch$site <- c("AZ", "ID", "KS", "MT", "NM")
corsynch <- round(cor(all_synch$mean_synch, all_synch$vravg),2)

R2 <- summary(lm(all_synch$mean_synch~all_synch$vravg))$r.squared
int <- round(coef(lm(all_synch$mean_synch~all_synch$vravg))[1],2)
slope <- round(coef(lm(all_synch$mean_synch~all_synch$vravg))[2],2)
# eq <- paste0("y = ", int, " + ", slope,"x")
eq <- expression(R^2)
library(ggthemes)
ggplot(all_synch, aes(x=vravg, y=mean_synch))+
  geom_abline(aes(intercept=0, slope=1), linetype=2)+
  geom_point()+
  geom_text(aes(label=site, y=mean_synch, x=vravg-0.04), size=3)+
  scale_x_continuous(limits=c(0,1))+
  scale_y_continuous(limits=c(0,1))+
  xlab("Average Synchrony of\nVital Rate Random Year Effects")+
  ylab("Average Synchrony of\nSimulated Per Capita Growth Rates")+
  theme_few()
ggsave(filename = "../docs/components/envresp_vs_ipmsynch.png", width = 4.5, height = 4, dpi = 150)







