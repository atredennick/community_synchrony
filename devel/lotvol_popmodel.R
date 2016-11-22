##  lotvol_popmodel.R: script to simulate population dynamics using
##  Lotka-Volterra model with environmental variability. Uses simulated
##  time series to calculate the coefficient of variation of total community
##  biomass with and without asynchronous environmental responses. Model assumes
##  interspecific competition is absent.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: November 14, 2016


rm(list=ls(all.names = TRUE))

set.seed(1234567)


####
####  Libraries ----------------------------------------------------------------
####
library(plyr)
library(reshape2)
library(mvtnorm)
library(synchrony)



####
####  Lotka-Volterra Model with Environmental Stochasticity --------------------
####
update_pop <- function(r, Nnow, K, env, sig_env){
  nspp <- length(Nnow)
  rm   <- numeric(nspp)
  for(i in 1:nspp){
    rm[i] <- (1 - (Nnow[i]/K[i])) + env[i]*sig_env[i]
  }
  return(rm)
}



####
####  Function to Generate Environmental Responses -----------------------------
####
get_env <- function(sigE, rho, nTime, num_spp) {
  varcov       <- matrix(rep(rho*sigE,num_spp*2), num_spp, num_spp)
  diag(varcov) <- sigE
  varcov <- as.matrix(varcov)
  e      <- rmvnorm(n = nTime, mean = rep(0,num_spp), sigma = varcov)
  return(e)
}



####
####  Simulate the Model -------------------------------------------------------
####
years_to_sim <- 2000
nspp <- 2
env_variance <- 1
rho <- seq(-1,1,by=0.05)
r <- rep(1, nspp)
K <- rep(1000, nspp)
sig_env <- rep(0.1, nspp)


cv_outs <- numeric(length(rho))
env_synch <- numeric(length(rho))

for(j in 1:length(rho)){
  fluct_env <- get_env(sigE = env_variance, rho = rho[j], nTime = years_to_sim, num_spp = nspp)
  N <- matrix(data = NA, nrow = years_to_sim, ncol = nspp)
  N[1,] <- 1
  rsaves <- matrix(data = NA, nrow = years_to_sim-1, ncol = nspp)
  
  for(t in 2:years_to_sim){
    rnows <- update_pop(r, N[t-1,], K, env = fluct_env[t,], sig_env)
    N[t,] <- N[t-1,] + N[t-1,]*rnows
    rsaves[t-1, ] <- rnows
  }
  
  matplot(N, type="l")
  cv <- sd(rowSums(N[500:2000,])) / mean(rowSums(N[500:2000,]))
  cv_outs[j] <- cv
  env_synch[j] <- as.numeric(community.sync(rsaves[500:1999,])[[1]])
}
plot(env_synch, cv_outs, frame.plot = FALSE, pch=19,
     xlab="Synchrony of Growth Rates", xlim=c(0,1),
     ylab="CV of Total Community Biomass")
cbind(env_synch, cv_outs)
# summary(lm(cv_outs~env_synch))
# abline(lm(cv_outs~env_synch), col="red")
# text(0.2,0.08,labels = paste("slope =",round(coef(lm(cv_outs~env_synch))[2],2)))
