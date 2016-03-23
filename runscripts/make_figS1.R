##  Main plot for synchrony results

library(ggplot2)
library(ggthemes)
library(reshape2)
library(plyr)
library(synchrony)
library(communitySynchrony)
library(gridExtra)

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

out_ks$percent_cover$site <- "Kansas"
out_id$percent_cover$site <- "Idaho"
out_mt$percent_cover$site <- "Montana"
out_nm$percent_cover$site <- "New Mexico"
out_az$percent_cover$site <- "Arizona"
cover_data <- rbind(out_ks$percent_cover, out_id$percent_cover,
                    out_mt$percent_cover, out_nm$percent_cover, 
                    out_az$percent_cover)

out_ks$growth_rates$site <- "Kansas"
out_id$growth_rates$site <- "Idaho"
out_mt$growth_rates$site <- "Montana"
out_nm$growth_rates$site <- "New Mexico"
out_az$growth_rates$site <- "Arizona"


site_colors <- c("grey45", "steelblue", "slateblue4", "darkorange", "purple")

nyrs <- c(max(out_az$stability$num_observations),
          max(out_id$stability$num_observations),
          max(out_ks$stability$num_observations),
          max(out_mt$stability$num_observations),
          max(out_nm$stability$num_observations))

boots <- 100
table_data <- readRDS("../results/calculate_site_metrics.RDS")

sites <- sort(table_data$site)
sites[5] <- "NewMexico"
synch <- matrix(NA, boots, length(sites))

output_list <- readRDS("../results/ipm_comp_nocomp_sims.RDS")
mlist <- melt(output_list)

for(do_site in 1:length(sites)){
  full_ipm_outs <- subset(mlist, L2=="ENVINTER")
  site_ipm <- subset(full_ipm_outs, L1==sites[do_site])
  colnames(site_ipm) <- c("year", "species", "cover", "L1", "L2")
  # if(sites[do_site]=="Idaho") site_ipm <- subset(site_ipm, Var2!="ARTR")
  species_list <- unique(site_ipm$species)
  num_spp <- length(species_list)
  for(i in 1:boots){
    begin_year <- sample(x = 1:(max(site_ipm$year)-nyrs[do_site]), 1)
    end_year <- begin_year+nyrs[do_site]
    tmp <- subset(site_ipm, year %in% begin_year:end_year)
    tmp2 <- dcast(tmp, year~species, value.var = "cover")
    
    # caclulate observed growth rates
    # create lagged data frame to only get observed yearly transitions
    lag_df <- tmp2
    lag_df$lagyear <- lag_df$year+1
    colnames(lag_df)[2:(num_spp+1)] <- paste0(colnames(lag_df)[2:(num_spp+1)],"_t0") 
    
    # merge the lag df with observed
    rm_col <- which(colnames(lag_df)=="year")
    merged_df <- merge(tmp2, lag_df[,-rm_col], by.x = "year", by.y="lagyear")
    transitions <- nrow(merged_df)
    obs_gr <- matrix(nrow=transitions, ncol=num_spp)
    name_ids1 <- which(colnames(merged_df) %in% species_list)
    name_ids2 <- which(colnames(merged_df) %in% paste(species_list, "_t0", sep=""))
    for(ii in 1:transitions){
      obs_gr[ii,] <- as.numeric(log(merged_df[ii,name_ids1]/merged_df[ii,name_ids2]))
    }
    synch[i,do_site] <- as.numeric(community.sync(obs_gr)[1])
  }
}

colnames(synch) <- sites
synch_df <- melt(synch)

obs_pgr_df <- c(as.numeric(out_az$pgr_synchrony[1]),
                as.numeric(out_id$pgr_synchrony[1]),
                as.numeric(out_ks$pgr_synchrony[1]),
                as.numeric(out_mt$pgr_synchrony[1]),
                as.numeric(out_nm$pgr_synchrony[1]))
obs_pgr_df <- data.frame(pgr_synch=obs_pgr_df,
                         Var2=sites)

### Get IBM synchs
synch_df$type <- "IPM"
synch_df <- synch_df[,2:ncol(synch_df)]
ibm_poly <- readRDS("../devel/mono_poly_ibm.RDS")
ibm_poly <- subset(ibm_poly, typesynch=="Per capita growth rate")
ibm_poly$type <- "IBM"
ibm_df <- ibm_poly[,c("site","synch","type")]
names(ibm_df) <- names(synch_df)
synch_df <- rbind(synch_df, ibm_df)

synch_df$site_order <- NA
synch_df[which(synch_df$Var2=="Arizona"), "site_order"] <- "2Arizona"
synch_df[which(synch_df$Var2=="Idaho"), "site_order"] <- "5Idaho"
synch_df[which(synch_df$Var2=="Kansas"), "site_order"] <- "3Kansas"
synch_df[which(synch_df$Var2=="Montana"), "site_order"] <- "4Montana"
synch_df[which(synch_df$Var2=="NewMexico"), "site_order"] <- "1NewMexico"

obs_pgr_df$site_order <- NA
obs_pgr_df[which(obs_pgr_df$Var2=="Arizona"), "site_order"] <- "2Arizona"
obs_pgr_df[which(obs_pgr_df$Var2=="Idaho"), "site_order"] <- "5Idaho"
obs_pgr_df[which(obs_pgr_df$Var2=="Kansas"), "site_order"] <- "3Kansas"
obs_pgr_df[which(obs_pgr_df$Var2=="Montana"), "site_order"] <- "4Montana"
obs_pgr_df[which(obs_pgr_df$Var2=="NewMexico"), "site_order"] <- "1NewMexico"



g1 <- ggplot()+
  geom_line(data=subset(synch_df, type=="IPM"), aes(x=value, color=site_order), 
            stat="density", adjust=2)+
  geom_vline(data=obs_pgr_df, aes(xintercept=pgr_synch, color=site_order), linetype=2)+
  facet_wrap("site_order", ncol=5, scales="free_y")+
  scale_color_manual(values = c("grey45", "steelblue", "slateblue4", "darkorange", "purple"))+
  guides(color=FALSE)+
  xlab("Community Synchrony")+
  ylab("Estimated Probability Density")+
  ggtitle("IPM")+
  scale_x_continuous(limits=c(0,1))+
  theme_bw()

g2 <- ggplot()+
  geom_line(data=subset(synch_df, type=="IBM"), aes(x=value, color=site_order), 
            stat="density", adjust=2)+
  geom_vline(data=obs_pgr_df, aes(xintercept=pgr_synch, color=site_order), linetype=2)+
  facet_wrap("site_order", ncol=5, scales="free_y")+
  scale_color_manual(values = c("grey45", "steelblue", "slateblue4", "darkorange", "purple"))+
  guides(color=FALSE)+
  xlab("Community Synchrony")+
  ylab("Estimated Probability Density")+
  ggtitle("IBM")+
  scale_x_continuous(limits=c(0,1))+
  theme_bw()

grid.arrange(g1,g2,nrow=2)


ggplot()+
  geom_line(data=subset(synch_df), aes(x=value, color=site_order), 
            stat="density", adjust=2)+
  geom_vline(data=obs_pgr_df, aes(xintercept=pgr_synch, color=site_order), linetype=2)+
  facet_grid(type~site_order, scales="free")+
  scale_color_manual(values = c("grey45", "steelblue", "slateblue4", "darkorange", "purple"))+
  guides(color=FALSE)+
  xlab("Community Synchrony")+
  ylab("Estimated Probability Density")+
  scale_x_continuous(limits=c(0,1))+
  theme_bw()

ggsave("../docs/components/figureS1.png", width = 14, height=4, units="in", dpi = 150)
