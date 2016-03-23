##
##  R Script to Simulate a Lotka-Volterra Model with Environmental Forcing
##
##  Author: Andrew Tredennick
##  Email: atredenn@gmail.com
##  Date created: 12-10-2015
##



####
####  Load Libraries
####
library(ggplot2) # for plotting
library(mvtnorm) # for multivariate normal


####
####  Set some global constants for model simulations
####
timesteps <- 1000      # number of iterations to run
r <- 0.25              # intrinsic growth rates
K <- 1000             # carrying capacities
Nstart <- K*2         # initial conditions
sigE <- 0.1           # environmental variance




####
####  Define Lotka-Volterra Model Function
####
####  Species are equivalent except for environmental response
####
update.lv <- function(n1, n2, r, K, evar1, evar2){
  r1 <- r * (1-(n1/K)) + evar1 
  r2 <- r * (1-(n2/K)) + evar2
  n1.new <- n1 + n1*r1
  n2.new <- n2 + n2*r2
  ntot <- n1.new+n2.new
  return(c(n1.new, n2.new, ntot))
}



####
####  Define function for (un)correlated environmental responses
####
get_env <- function(sigE, rho_spp, timesteps){
  varcov <- matrix(c(sigE, rho_spp*sigE,
                     rho_spp*sigE, sigE), 
                   2, 2)
  e <- rmvnorm(n = timesteps, mean = c(0,0), sigma = varcov)
  # g <- exp(e) / (1+exp(e))
  return(e)
}



####
####  Simulate Model
####
rho.vec <- seq(-1,1,length.out = 100)
cv.out <- numeric(length(rho.vec))
e.synch <- numeric(length(rho.vec))
for(i in 1:length(rho.vec)){
  rho_spp <- rho.vec[i]
  evar <- get_env(1, rho_spp, timesteps)
  evar <- evar*sigE
  n.pop <- matrix(NA, nrow=timesteps, ncol=3)
  n.pop[1,1:2] <- Nstart
  n.pop[1,3] <- Nstart*2
  for(t in 2:timesteps){
    n.pop[t,] <- update.lv(n.pop[t-1,1], n.pop[t-1,2], r, K, evar[t,1], evar[t,2])
  }
  
  # matplot(n.pop, type="l")
  cv.out[i] <- sd(n.pop[20:timesteps,3])/mean(n.pop[20:timesteps,3])
  e.synch[i] <- as.numeric(community.sync(evar)[1])
}




plot_df <- data.frame(rho=rho.vec, cv=cv.out)
ggplot(plot_df, aes(x=rho, y=cv))+
  geom_point(size=3, shape=1)+
  stat_smooth(se=FALSE, size=1, color="black", method="loess")+
  # stat_smooth(se=FALSE, size=1.2, color="black", linetype=2, method="lm")+
  xlab(bquote(rho[e]))+
  ylab(bquote(CV[NT]))+
  theme_bw()
ggsave("../results/lotka_volt_cv_rho.png", width = 3.5, height=3)
# geom_point(size=3, color="#9D6188", alpha=0.7)+
#   stat_smooth(se=FALSE, color="#97A761", size=1.5)+
