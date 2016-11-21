##
##  Development Script to Look at Individual Simulated Time Series
##

rm(list=ls())

# par(mfrow=c(1,2))
ibms <- readRDS("../results/ibm_sims/ibm_Kansas_constnointerexpand1.RDS")
matplot(ibms[,c(4:5)]*100, type="l")


mycols <- c("#BACB72", "#1E4B58", "#3997B3", "#57A880")

ibms <- readRDS("../results/ibm_sims/ibm_Montana_fluctnointerexpand5.RDS")
matplot(ibms[101:200,c(4:7)]*100, type="l", lty = 1, lwd=2, col = mycols, bty="n", xaxt='n', yaxt='n', ann=FALSE)
axis(1, labels = FALSE)
axis(2, labels = FALSE)

# caclulate observed growth rates
# create lagged data frame to only get observed yearly transitions
ts_mat <- as.data.frame(ibms[1:100,c(2,4:7)])
lag_df <- ts_mat
lag_df$lagyear <- lag_df$time+1
num_spp=4
species_list <- c("Cov.BOGR","Cov.HECO","Cov.PASM","Cov.POSE")
colnames(lag_df)[2:(num_spp+1)] <- paste(colnames(lag_df)[2:(num_spp+1)],"_t0", sep="") 
# merge the lag df with observed
rm_col <- which(colnames(lag_df)=="time")
merged_df <- merge(ts_mat, lag_df[,-rm_col], by.x = "time", by.y="lagyear")
transitions <- nrow(merged_df)
obs_gr <- matrix(nrow=transitions, ncol=num_spp)
name_ids1 <- which(colnames(merged_df) %in% species_list)
name_ids2 <- which(colnames(merged_df) %in% paste(species_list, "_t0", sep=""))
for(i in 1:transitions){
  obs_gr[i,] <- as.numeric(log(merged_df[i,name_ids1]/merged_df[i,name_ids2]))
}
matplot(obs_gr, type="l", lty = 1, lwd=2, col = mycols, bty="n", xaxt='n', yaxt='n', ann=FALSE)
axis(1, labels = FALSE)
axis(2, labels = FALSE)
