##
##  Development Script to Look at Individual Simulated Time Series
##

rm(list=ls())

par(mfrow=c(1,2))
ibms <- readRDS("../results/ibm_sims/ibm_Kansas_constnointerexpand1.RDS")
matplot(ibms[,c(4:5)]*100, type="l")
ibms <- readRDS("../results/ibm_sims/ibm_Kansas_constnointerexpand5.RDS")
matplot(ibms[,c(4:5)]*100, type="l")
