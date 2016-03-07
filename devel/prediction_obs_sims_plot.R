##  R script to plot comparisons between observed/simulated synchrony and
##  expectations under two limiting conditions

rm(list=ls())

library(ggplot2)
library(ggthemes)
library(plyr)
library(reshape2)
library(communitySynchrony)
library(synchrony)

####
#### Observed synchrony (cover and growth rates)
####
site <- "Kansas"
spp_list <- c("BOCU","BOHI","SCSC")
num_spp <- length(spp_list)
ks_data <- data.frame(quad=NA, year=NA, totCover=NA, species=NA)
for(dospp in 1:num_spp){ #loop through species to read in data
  spp_now <- spp_list[dospp]
  quad_file <- paste("../data/", site,"/",spp_now,"/quadratCover.csv",sep="")
  spp_data <- read.csv(quad_file)
  spp_data$species <- spp_now
  ks_data <- rbind(ks_data, spp_data)
} #end species looping for raw data
ks_data <- ks_data[2:nrow(ks_data),] #remove first NA row

tmp1<-which(ks_data$quad_data=="q25" & (ks_data$year<35 | ks_data$year>62))
tmp2<-which(ks_data$quad_data=="q27")
tmp3<-which(ks_data$quad=="q28")
tmp4<-which(ks_data$quad=="q30")
tmp5<-which(ks_data$quad=="q31" & (ks_data$year<35 | ks_data$year>39))
tmp6<-which(ks_data$quad=="q32" & (ks_data$year<35 | ks_data$year>41))
tmp<-c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
ks_data<-ks_data[-tmp,]

# exclude the records later than 1968, to keep the same random year effect...
ks_data<-subset(ks_data,year<68)


site <- "Idaho"
spp_list <- sort(c("PSSP","HECO","POSE","ARTR"))
num_spp <- length(spp_list)
id_data <- data.frame(quad=NA, year=NA, totCover=NA, species=NA)
for(dospp in 1:num_spp){ #loop through species to read in data
  spp_now <- spp_list[dospp]
  quad_file <- paste("../data/", site,"/",spp_now,"/quadratCover.csv",sep="")
  spp_data <- read.csv(quad_file)
  spp_data$species <- spp_now
  id_data <- rbind(id_data, spp_data)
} #end species looping for raw data
id_data <- id_data[2:nrow(id_data),] #remove first NA row
# id_data <- subset(id_data, species!="ARTR") #take out the shrub


site <- "Montana"
spp_list <- sort(c("BOGR","HECO","PASM","POSE"))
num_spp <- length(spp_list)
mt_data <- data.frame(quad=NA, year=NA, totCover=NA, species=NA)
for(dospp in 1:num_spp){ #loop through species to read in data
  spp_now <- spp_list[dospp]
  quad_file <- paste("../data/", site,"/",spp_now,"/quadratCover.csv",sep="")
  spp_data <- read.csv(quad_file)
  spp_data$species <- spp_now
  mt_data <- rbind(mt_data, spp_data)
} #end species looping for raw data
mt_data <- mt_data[2:nrow(mt_data),] #remove first NA row



site <- "NewMexico"
spp_list <- sort(c("BOER","SPFL"))
num_spp <- length(spp_list)
nm_data <- data.frame(quad=NA, year=NA, totCover=NA, species=NA)
for(dospp in 1:num_spp){ #loop through species to read in data
  spp_now <- spp_list[dospp]
  quad_file <- paste("../data/", site,"/",spp_now,"/quadratCover.csv",sep="")
  spp_data <- read.csv(quad_file)
  spp_data$species <- spp_now
  nm_data <- rbind(nm_data, spp_data)
} #end species looping for raw data
nm_data <- nm_data[2:nrow(nm_data),] #remove first NA row



site <- "Arizona"
spp_list <- sort(c("BOER","BORO"))
num_spp <- length(spp_list)
az_data <- data.frame(quad=NA, year=NA, totCover=NA, species=NA)
for(dospp in 1:num_spp){ #loop through species to read in data
  spp_now <- spp_list[dospp]
  quad_file <- paste("../data/", site,"/",spp_now,"/quadratCover.csv",sep="")
  spp_data <- read.csv(quad_file)
  spp_data$species <- spp_now
  az_data <- rbind(az_data, spp_data)
} #end species looping for raw data
az_data <- az_data[2:nrow(az_data),] #remove first NA row

out_ks <- get_comm_synchrony(ts_data = ks_data)
out_id <- get_comm_synchrony(ts_data = id_data)
out_mt <- get_comm_synchrony(ts_data = mt_data)
out_nm <- get_comm_synchrony(ts_data = nm_data)
out_az <- get_comm_synchrony(ts_data = az_data)

get_table_metrics <- function(output){
  tmp <- output
  exp_tmp <- round(tmp$pgr_expected_synch_ind_flucts, 2)
  obs_tmp <- round(as.numeric(tmp$pgr_synchrony[1]), 2)
  demo_tmp <- obs_tmp - exp_tmp
  exp_cover <- round(as.numeric(tmp$cover_expected_synch_ind_flucts), 2)
  obs_cover <- round(as.numeric(tmp$abund_synchrony[1]), 2)
  demo_cover <- obs_cover - exp_cover
  plant_size <- tmp$avg_size
  cover_var <- tmp$variability
  return(c(obs_tmp, exp_tmp, demo_tmp, obs_cover, exp_cover, 
           demo_cover, cover_var))
}

table_data <- as.data.frame(rbind(get_table_metrics(out_az),
                                  get_table_metrics(out_ks),
                                  get_table_metrics(out_id),
                                  get_table_metrics(out_mt),
                                  get_table_metrics(out_nm)))

site_names <- c("Arizona", "Kansas", "Idaho", "Montana", "NewMexico")
obs_synchronies <- table_data[,c("V1", "V4")]
obs_synchronies$site <- site_names
colnames(obs_synchronies) <- c("pgr", "cover", "site")
exp_synchrony_demo <- table_data[,c("V2","V5")]
exp_synchrony_demo$site <- site_names
colnames(exp_synchrony_demo) <- c("pgr", "cover", "site")


####
#### Simulated synchrony (cover and growth rates)
####
output_list <- readRDS("../results/ipm_comp_nocomp_sims.RDS")
mlist <- melt(output_list)
colnames(mlist)[1:3] <- c("year", "species", "cover")
sites <- unique(mlist$L1)
sims <- unique(mlist$L2)
synch_df <- data.frame(site=NA, experiment=NA, bootnum=NA, 
                       pgr_synch=NA, abund_synch=NA)
boots <- 100
num_iters <- 50
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
agg_synch <- subset(agg_synch, experiment=="ENVINTER")



####
#### Expected synchrony (cover and growth rates; two limiting cases)
####
exp_synchrony_env <- readRDS("../results/envstoch_predictions.RDS")
exp_synchrony_env <- exp_synchrony_env[complete.cases(exp_synchrony_env),]


####
####  Put them all together and plot
####
plot_df <- data.frame(site=c(rep(site_names,4)),
                      obs_synchrony=c(obs_synchronies$cover,
                                      obs_synchronies$pgr,
                                      obs_synchronies$cover,
                                      obs_synchronies$pgr),
                      pred_synchrony=c(exp_synchrony_env$cover_prediction,
                                       exp_synchrony_env$growthrate_prediction,
                                       exp_synchrony_demo$cover,
                                       exp_synchrony_demo$cover),
                      typesynch=c(rep("Percent Cover", 5),
                                  rep("Growth Rate",5),
                                  rep("Percent Cover", 5),
                                  rep("Growth Rate",5)),
                      predtype=c(rep("env",5*2),
                                 rep("demo",5*2)))


mycolors <- c("#9D6188","#97A861")
gout <- ggplot(subset(plot_df, typesynch=="Growth Rate"), aes(x=pred_synchrony, y=obs_synchrony, color=predtype))+
  geom_point(size=3)+
  stat_smooth(method="lm", se=FALSE, size=1)+
  geom_abline(aes(intercept=0, slope=1), linetype=2, color="grey")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(0,1))+
  # facet_wrap("typesynch", ncol=1)+
  scale_color_manual(values=mycolors, labels=c(expression(M[D]), expression(M[E])), name="")+
  xlab("Predicted Synchrony")+
  ylab("Observed Synchrony")+
  # guides(shape=FALSE)+
  theme_few()+
  theme(legend.position=c(0.1,0.85),
        legend.background = element_rect(colour = NA, fill = NA))

png("../docs/components/prediction_observed.png", width = 4, height=4, units = "in", res=100)
print(gout)
dev.off()

corrs_all <- ddply(plot_df, .(typesynch, predtype), summarise,
                   value = cor(pred_synchrony, obs_synchrony))

env_pgr_preds <- predict(lm(obs_synchrony~pred_synchrony, data=subset(plot_df, typesynch=="growth rate" & predtype=="env")))
demo_pgr_preds <- predict(lm(obs_synchrony~pred_synchrony, data=subset(plot_df, typesynch=="growth rate" & predtype=="demo")))
env_cover_preds <- predict(lm(obs_synchrony~pred_synchrony, data=subset(plot_df, typesynch=="cover" & predtype=="env")))
demo_cover_preds <- predict(lm(obs_synchrony~pred_synchrony, data=subset(plot_df, typesynch=="cover" & predtype=="demo")))

errors_env_pgr <- unlist(subset(plot_df, typesynch=="growth rate" & predtype=="env")["pred_synchrony"])-unlist(subset(plot_df, typesynch=="growth rate" & predtype=="env")["obs_synchrony"])
preds <- unlist(subset(plot_df, typesynch=="growth rate" & predtype=="env")["pred_synchrony"])
summary(lm(errors_env_pgr~preds))
errors_dem_pgr <- unlist(subset(plot_df, typesynch=="growth rate" & predtype=="demo")["pred_synchrony"])-unlist(subset(plot_df, typesynch=="growth rate" & predtype=="demo")["obs_synchrony"])
preds <- unlist(subset(plot_df, typesynch=="growth rate" & predtype=="demo")["pred_synchrony"])
summary(lm(errors_dem_pgr~demo_pgr_preds))
