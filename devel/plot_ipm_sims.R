library(ggplot2)
library(reshape2)
library(plyr)
library(synchrony)
library(communitySynchrony)

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

ipm_pgr_plot <- ggplot(agg_synch, aes(x=experiment, y=mean_synch, 
                                      color=site, group=typesynch))+
  geom_line(aes(linetype=typesynch))+
  geom_errorbar(aes(ymin=lo_synch, ymax=up_synch, linetype=typesynch), width=0.05)+
  geom_point(size=2, aes(shape=typesynch))+
  scale_color_manual(values = c("grey45", "steelblue", "slateblue4", "darkorange", "purple"),
                     name = "Site")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_discrete(labels=c("F-INTER", "F-NoINTER"))+
  scale_linetype_discrete(name="", labels=c("Per capita growth rate synchrony",
                                            "Cover (%) synchrony"))+
  scale_shape_discrete(name="", labels=c("Per capita growth rate synchrony",
                                            "Cover (%) synchrony"))+
  xlab("Simulation")+
  ylab("Community synchrony")+
  facet_wrap("site", nrow=1)+
  guides(color=FALSE)+
  theme_bw()+
  theme(legend.position = "top",
        legend.background = element_rect(fill=NA,
                                         size=0.5))

png("../docs/components/ipm_pgr_plot.png", width=8.5, height=3, units="in", res=150)
print(ipm_pgr_plot)
dev.off()



####
####  Plot synchrony in polyculture vs synchrony in monoculture
####
library(ggthemes)
polymono <- agg_synch[,c("site", "typesynch", "experiment", "mean_synch")]
polymono_wide <- dcast(polymono, site+typesynch~experiment, value.var = "mean_synch")
g1 <- ggplot(polymono_wide, aes(x=ENVNOINTER, y=ENVINTER))+
  geom_abline(aes(intercept=0, slope=1), linetype=3)+
  geom_point(size=3, aes(shape=typesynch, color=site))+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(0,1))+
  ylab("Species synchrony in polyculture")+
  xlab("Species synchrony \nin monoculture")+
  scale_shape_discrete(name="Temporal Variable", labels=c("Per capita growth rate", "Percent cover"))+
  scale_color_manual(values = c("grey45", "steelblue", "slateblue4", "darkorange", "purple"),
                     name = "Site")+
  theme_few()+
  guides(shape=FALSE)+
  ggtitle("A                                                          ")+
  theme(legend.position=c(0.2,0.7))+
  theme(legend.text = element_text(size = 8))+
  theme(legend.background = element_rect(colour = "grey", fill = NA))



output_list <- readRDS("../results/ipm_yearly_pgr.RDS")
mlist <- melt(output_list)
colnames(mlist)[1:3] <- c("year", "species", "pgr")
sites <- unique(mlist$L1)
synch_df <- data.frame(site=NA, bootnum=NA, pgr_synch=NA)
boots <- 100
num_iters <- 50
for(dosite in sites){
  tmpsim <- subset(mlist, L1==dosite)
    for(i in 1:boots){
      begin_year <- sample(x = 1:(max(tmpsim$year)-num_iters), 1)
      end_year <- begin_year+num_iters
      tmp <- subset(tmpsim, year %in% begin_year:end_year)
      tmp <- dcast(tmp, year~species, value.var = "pgr")
      tmpsynch <- as.numeric(community.sync(tmp[2:ncol(tmp)])[1])
      tmp_df <- data.frame(site=dosite, bootnum=i, pgr_synch=tmpsynch)
      synch_df <- rbind(synch_df, tmp_df)
    }# end boots loop
}# end site loop
synch_dftmp <- synch_df[2:nrow(synch_df),]
pgr_synch <- ddply(synch_dftmp, .(site), summarise,
                   mean_pgrsynch = mean(pgr_synch))



# pgr_list <- readRDS("../results/ipm_yearly_pgr.RDS")
# site_names <- names(pgr_list)
# out_df <- data.frame(site=NA, pgr_synch=NA)
# for(do_site in site_names){
#   tmp_pgrs <- pgr_list[[do_site]]
#   tmp_synch <- as.numeric(community.sync(tmp_pgrs)[1])
#   tmp_df <- data.frame(site=do_site, pgr_synch=tmp_synch)
#   out_df <- rbind(out_df, tmp_df)
# }
# pgr_synch <- out_df[2:nrow(out_df),]
polymono_wide$yrpgr <- rep(pgr_synch$mean_pgrsynch, each=2)

g2 <- ggplot(polymono_wide, aes(x=yrpgr, y=ENVINTER))+
  geom_abline(aes(intercept=0, slope=1), linetype=3)+
  geom_point(size=3, aes(shape=typesynch, color=site))+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(0,1))+
  ylab("Species synchrony in polyculture")+
  xlab("Yearly per capita growth rate synchrony \nin monoculture")+
  scale_shape_discrete(name="Temporal Variable", labels=c("Per capita growth rate", "Percent cover"))+
  scale_color_manual(values = c("grey45", "steelblue", "slateblue4", "darkorange", "purple"),
                     name = "Site")+
  theme_few()+
  guides(shape=FALSE, color=FALSE)+
  ggtitle("B                                                          ")

library(gridExtra)

png("../docs/components/poly_vs_mono_synch.png",width = 8, height = 4, units = "in", res=150)
g_out <- grid.arrange(g1,g2,ncol=2)
dev.off()

# ggsave("../docs/components/poly_vs_mono_synch.png", plot=g_out, width = 8, height = 3.5, dpi = 150)
