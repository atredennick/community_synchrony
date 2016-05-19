##  Main plot for synchrony results

library(ggplot2)
library(ggthemes)
library(reshape2)
library(plyr)
library(synchrony)
library(communitySynchrony)

site_colors <- c("grey45", "steelblue", "slateblue4", "darkorange", "purple")



####
####  MAIN TEXT FIGURES 1-2
####
##  Read in IPM results ---------
output_list <- readRDS("../results/ipm_comp_nocomp_sims.RDS")
mlist <- melt(output_list)
colnames(mlist)[1:3] <- c("year", "species", "cover")
sites <- unique(mlist$L1)
sims <- unique(mlist$L2)
synch_df <- list()

num_iters <- 50
site_species_variances <- list()
png("../docs/components/environmental_variances.png", width = 10, height=3, units="in", res=100)
par(mfrow=c(1,5), las=3)
for(dosite in sites){
  tmp_data <- subset(mlist, L1==dosite)
  for(dosim in sims){
    tmpsim <- subset(tmp_data, L2==dosim)
    tmp_abund_vars <- matrix(ncol=length(unique(tmpsim$species)), nrow=boots)
    for(i in 1:boots){
      begin_year <- sample(x = 1:(max(tmpsim$year)-num_iters), 1)
      end_year <- begin_year+num_iters
      tmp <- subset(tmpsim, year %in% begin_year:end_year)
      tmpcast <- dcast(tmp, year~species, value.var = "cover")
      tmp_abund_vars[i,] <- (apply(tmpcast[2:ncol(tmpcast)], 2, "sd"))^2
    }# end boots loop
    site_species_variances[[dosite]][[dosim]] <- tmp_abund_vars
    colnames(tmp_abund_vars) <- unique(tmpsim$species)
    if(dosim=="ENVNOINTER") { boxplot(tmp_abund_vars, main=dosite,  ylab="environmental variance", outline=F) }
  }# end experiment/sim loop
}# end site loop
dev.off()







