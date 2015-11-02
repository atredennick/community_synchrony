library(ggplot2)
library(reshape2)
library(plyr)
library(communitySynchrony)

output_list <- readRDS("../results/ipm_comp_nocomp_sims.RDS")
mlist <- melt(output_list)
colnames(mlist)[1:3] <- c("year", "species", "cover")
sites <- unique(mlist$L1)
sims <- unique(mlist$L2)
synch_df <- data.frame(site=NA, experiment=NA, bootnum=NA, synch=NA)
boots <- 100
num_iters <- 100
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
      tmp_df <- data.frame(site=dosite, experiment=dosim, 
                           bootnum=i, synch=tmp_pgr_synch)
      synch_df <- rbind(synch_df, tmp_df)
    }# end boots loop
  }# end experiment/sim loop
}# end site loop

synch_df <- synch_df[2:nrow(synch_df),]
agg_synch <- ddply(synch_df, .(site, experiment), summarise,
                   mean_synch = mean(synch),
                   up_synch = quantile(synch, 0.95),
                   lo_synch = quantile(synch, 0.05))

ipm_pgr_plot <- ggplot(agg_synch, aes(x=experiment, y=mean_synch, color=site, group=site))+
  geom_line()+
  geom_errorbar(aes(ymin=lo_synch, ymax=up_synch), width=0.05)+
  geom_point(size=2)+
  scale_color_manual(values = c("grey45", "steelblue", "slateblue4", "darkorange", "purple"),
                     name = "Site")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_discrete(labels=c("F-INTER", "F-NoINTER"))+
  xlab("Simulation")+
  ylab("Species synchrony")+
  facet_wrap("site", nrow=1)+
  guides(color=FALSE)+
  theme_bw()

png("../docs/components/ipm_pgr_plot.png", width=8.5, height=2.5, units="in", res=200)
print(ipm_pgr_plot)
dev.off()
