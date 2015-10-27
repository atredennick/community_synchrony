library(communitySynchrony)

output_list <- readRDS("../results/ipm_comp_nocomp_sims.RDS")
mlist <- melt(output_list)
colnames(mlist)[1:3] <- c("year", "species", "cover")
sites <- unique(mlist$L1)
sims <- unique(mlist$L2)
pgr_list <- list()
for(dosite in sites){
  tmp_data <- subset(mlist, L1==dosite)
  if(dosite=="Idaho"){
    tmp_data <- subset(tmp_data, species!="ARTR")
  }
  within_site_list <- list()
  for(dosim in sims){
    tmpsim <- subset(tmp_data, L2==dosim)[c("year", "species", "cover")]
    tmpsynch <- get_ipm_synchrony(tmpsim)
    tmp_pgr_synch <- as.numeric(tmpsynch$pgr_synchrony["obs"])
    within_site_list[[dosim]] <- tmp_pgr_synch
  }
  pgr_list[[dosite]] <- within_site_list
}
pgr_synch <- melt(pgr_list)
expected_synch <- table_data[,c(1,3)]
expected_synch$L2 <- "E(SYNCH)"
colnames(expected_synch) <- c("L1", "value", "L2")
expected_synch <- expected_synch[,c("value", "L2", "L1")]
expected_synch[which(expected_synch[,"L1"]=="New Mexico"),"L1"] <- "NewMexico"

pgr_synch <- rbind(pgr_synch, expected_synch)

# ggplot(pgr_synch, aes(x=L2, y=value, fill=L2))+
#   geom_bar(stat="identity", position=position_dodge(0.9), color="grey")+
#   facet_wrap("L1", scales="free_x")+
#   scale_fill_manual(values = c("grey45", "steelblue", "slateblue4"))+
#   guides(fill=FALSE)+
#   scale_y_continuous(limits=c(0,1))+
#   xlab("Source of stochasticity")+
#   ylab("Species synchrony")

ggplot(pgr_synch, aes(x=L2, y=value, color=L1, group=L1))+
  geom_line()+
  geom_point(size=5)+
  scale_color_manual(values = c("grey45", "steelblue", "slateblue4", "darkorange", "purple"),
                     name = "Site")+
  scale_y_continuous(limits=c(0,1))+
  xlab("Simulation")+
  ylab("Species synchrony")+
  theme_bw()