##  Functions using INLA to get vital rate regression parameters

#' Estimate growth regression coefficients using INLA
#' 
#' @author Andrew Tredennick
#' @param dataframe Time series dataframe of genet sizes at t0 and t1 with 
#'                  neighborhood crowding covariate for each observation.
#' @param crowd_mat Matrix of crowding indices for each observation. Matrix
#'                  dimensions are nrow = number of observations,
#'                  ncol = number of species.
#' @param alphas Vector of length(species) of vital rate specific alpha values.
#' @return Dataframe with named regression coefficients.

get_growth_params_yrlcomp <- function(dataframe, crowd_mat, alpha){
  D <- dataframe
  D$logarea.t0 <- log(D$area.t0)
  D$logarea.t1 <- log(D$area.t1)
  D$quad <- as.character(D$quad)
  
  crowd = crowd_mat 
  crowd[crowd<1e-99]=0 # make really small crowding indices 0
  
  library(INLA)
  # Set up ID variables for INLA random effects
  D$yearID <- D$year+max(D$year) # for random year offset on intercept
  D$GroupID <- as.numeric(D$Group)
  
  if(ncol(crowd)==2) {
    W1 <- crowd[,1]
    W2 <- crowd[,2]
    D$yearW1 <- 100+D$yearID # for random year effect on crowding
    D$yearW2 <- 100+D$yearW1 # for random year effect on crowding
    
    formula <- logarea.t1 ~ logarea.t0+W1+W2+
      f(yearID, model="iid", prior="normal",param=c(0,0.001))+
      f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
      f(year, logarea.t0, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW1, W1, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW2, W2, model="iid", prior="normal",param=c(0,0.001))
  }
  if(ncol(crowd)==3) {
    W1 <- crowd[,1]
    W2 <- crowd[,2]
    W3 <- crowd[,3]
    D$yearW1 <- 100+D$yearID # for random year effect on crowding
    D$yearW2 <- 100+D$yearW1 # for random year effect on crowding
    D$yearW3 <- 100+D$yearW2 # for random year effect on crowding
    
    formula <- logarea.t1 ~ logarea.t0+W1+W2+W3+
      f(yearID, model="iid", prior="normal",param=c(0,0.001))+
      f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
      f(year, logarea.t0, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW1, W1, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW2, W2, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW3, W3, model="iid", prior="normal",param=c(0,0.001))
  }
  if(ncol(crowd)==4) {
    W1 <- crowd[,1]
    W2 <- crowd[,2]
    W3 <- crowd[,3]
    W4 <- crowd[,4]
    D$yearW1 <- 100+D$yearID # for random year effect on crowding
    D$yearW2 <- 100+D$yearW1 # for random year effect on crowding
    D$yearW3 <- 100+D$yearW2 # for random year effect on crowding
    D$yearW4 <- 100+D$yearW3 # for random year effect on crowding
    
    formula <- logarea.t1 ~ logarea.t0+W1+W2+W3+W4+
      f(yearID, model="iid", prior="normal",param=c(0,0.001))+
      f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
      f(year, logarea.t0, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW1, W1, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW2, W2, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW3, W3, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW4, W4, model="iid", prior="normal",param=c(0,0.001))
  }
  
  outINLA <- inla(formula, data=D,
                  family=c("gaussian"), verbose=FALSE,
                  control.predictor = list(link = 1),
                  control.compute=list(dic=T,mlik=T),
                  control.inla = list(h = 1e-10),
                  control.fixed = list(correlation.matrix=TRUE))
  rand_eff_summaries <- outINLA$summary.random
#   tmp_fixed_covariance <- outINLA$misc$lincomb.derived.correlation.matrix
#   
#   # Fit variance
#   x <- outINLA$summary.fitted.values$mean #fitted values from INLA
#   y <- (D$logarea.t1-outINLA$summary.fitted.values$mean)^2 #calculates the variance (residuals^2)
#   outVar <- nls(y~a*exp(b*x),start=list(a=1,b=0)) #model the variance as function of fitted size
  
  # Collect parameters
  # random year and group effects
  if(ncol(crowd)==2){
    params <- as.data.frame(outINLA$summary.random$yearID[,1:2])
    params <- cbind(params, outINLA$summary.random$year[,2])
    params <- cbind(params, outINLA$summary.random$yearW1[,2])
    params <- cbind(params, outINLA$summary.random$yearW2[,2])
    names(params) <- c("Year", "Intercept.yr", "logarea.t0.yr", 
                       "W1.yr", "W2.yr")
  }
  if(ncol(crowd)==3){
    params <- as.data.frame(outINLA$summary.random$yearID[,1:2])
    params <- cbind(params, outINLA$summary.random$year[,2])
    params <- cbind(params, outINLA$summary.random$yearW1[,2])
    params <- cbind(params, outINLA$summary.random$yearW2[,2])
    params <- cbind(params, outINLA$summary.random$yearW3[,2])
    names(params) <- c("Year", "Intercept.yr", "logarea.t0.yr", 
                       "W1.yr", "W2.yr", "W3.yr")
  }
  if(ncol(crowd)==4){
    params <- as.data.frame(outINLA$summary.random$yearID[,1:2])
    params <- cbind(params, outINLA$summary.random$year[,2])
    params <- cbind(params, outINLA$summary.random$yearW1[,2])
    params <- cbind(params, outINLA$summary.random$yearW2[,2])
    params <- cbind(params, outINLA$summary.random$yearW3[,2])
    params <- cbind(params, outINLA$summary.random$yearW4[,2])
    names(params) <- c("Year", "Intercept.yr", "logarea.t0.yr", 
                       "W1.yr", "W2.yr", "W3.yr", "W4.yr")
  }
  
  tmp <- as.data.frame(outINLA$summary.random$GroupID[,2])
  names(tmp) <- "Group"
  tmp$GrpName <- outINLA$summary.random$GroupID[,1]
  tmp[(NROW(tmp)+1):NROW(params),]=NA
  params=cbind(params,tmp)
  
  # fixed effects
  fixed <- as.data.frame(outINLA$summary.fixed)[,1:2]
  tmp=matrix(NA,dim(params)[1],nrow(fixed))
  colnames(tmp)=rownames(fixed)
  tmp[1,]=fixed[,1]
  params=cbind(params,tmp)
  params$alpha=NA; params$alpha[1:length(alpha)]=alpha
  colnames(params)[which(colnames(params)=="(Intercept)")] <- "Intercept"
#   #variance 
#   params$sigma.a=NA; params$sigma.a[1]=coef(outVar)[1] 
#   params$sigma.b=NA; params$sigma.b[1]=coef(outVar)[2]
  
  return(list(params=params, rand_eff_summaries=rand_eff_summaries))
}






#' Estimate survival regression coefficients using INLA
#' 
#' @author Andrew Tredennick
#' @param dataframe Time series dataframe of genet sizes at t0 and t1 with 
#'                  neighborhood crowding covariate for each observation.
#' @param crowd_mat Matrix of crowding indices for each observation. Matrix
#'                  dimensions are nrow = number of observations,
#'                  ncol = number of species.
#' @param alphas Vector of length(species) of vital rate specific alpha values.
#' @return Dataframe with named regression coefficients.

get_survival_params_yrlcomp <- function(dataframe, crowd_mat, alpha){
  D <- dataframe
  D$logarea <- log(D$area)
  D$quad <- as.character(D$quad)
  
  crowd = crowd_mat 
  crowd[crowd<1e-99]=0 # make really small crowding indices 0
  
  library(INLA)
  # Set up ID variables for INLA random effects
  D$yearID <- D$year+max(D$year) # for random year offset on intercept
  D$GroupID <- as.numeric(D$Group)
  
  if(ncol(crowd)==2){
    W1 <- crowd[,1]
    W2 <- crowd[,2]
    D$yearW1 <- 100+D$yearID # for random year effect on crowding
    D$yearW2 <- 100+D$yearW1 # for random year effect on crowding
    
    formula2 <- survives ~ logarea+W1+W2+
      f(yearID, model="iid", prior="normal",param=c(1,0.001))+
      f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
      f(year, logarea, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW1, W1, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW2, W2, model="iid", prior="normal",param=c(0,0.001))
  }
  if(ncol(crowd)==3){
    W1 <- crowd[,1]
    W2 <- crowd[,2]
    W3 <- crowd[,3]
    D$yearW1 <- 100+D$yearID # for random year effect on crowding
    D$yearW2 <- 100+D$yearW1 # for random year effect on crowding
    D$yearW3 <- 100+D$yearW2 # for random year effect on crowding
    
    formula2 <- survives ~ logarea+W1+W2+W3+
      f(yearID, model="iid", prior="normal",param=c(1,0.001))+
      f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
      f(year, logarea, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW1, W1, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW2, W2, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW3, W3, model="iid", prior="normal",param=c(0,0.001))
  }
  if(ncol(crowd)==4){
    W1 <- crowd[,1]
    W2 <- crowd[,2]
    W3 <- crowd[,3]
    W4 <- crowd[,4]
    D$yearW1 <- 100+D$yearID # for random year effect on crowding
    D$yearW2 <- 100+D$yearW1 # for random year effect on crowding
    D$yearW3 <- 100+D$yearW2 # for random year effect on crowding
    D$yearW4 <- 100+D$yearW3 # for random year effect on crowding
    
    formula2 <- survives ~ logarea+W1+W2+W3+W4+
      f(yearID, model="iid", prior="normal",param=c(1,0.001))+
      f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
      f(year, logarea, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW1, W1, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW2, W2, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW3, W3, model="iid", prior="normal",param=c(0,0.001))+
      f(yearW4, W4, model="iid", prior="normal",param=c(0,0.001))
  }
  
  outINLA <- inla(formula2, data=D,
                  family=c("binomial"), verbose=FALSE,
                  control.compute=list(dic=T,mlik=T),
                  control.predictor = list(link = 1),
                  control.inla = list(h = 1e-10),
                  Ntrials=rep(1,nrow(D)),
                  control.fixed = list(correlation.matrix=TRUE))
  rand_eff_summaries <- outINLA$summary.random
  # tmp_fixed_covariance <- outINLA$misc$lincomb.derived.correlation.matrix
  
  #Collect parameters
  #random year and group effects
  if(ncol(crowd)==2){
    params <- as.data.frame(outINLA$summary.random$yearID[,1:2])
    params <- cbind(params, outINLA$summary.random$year[,2])
    params <- cbind(params, outINLA$summary.random$yearW1[,2])
    params <- cbind(params, outINLA$summary.random$yearW2[,2])
    ngroups <- length(outINLA$summary.random$GroupID[,2])
    params <- cbind(params, c(outINLA$summary.random$GroupID[,2],rep(NA,nrow(params)-ngroups)))
    names(params) <- c("Year", "Intercept.yr", "logarea.t0.yr", 
                       "W1.yr", "W2.yr", "Group")
  }
  if(ncol(crowd)==3){
    params <- as.data.frame(outINLA$summary.random$yearID[,1:2])
    params <- cbind(params, outINLA$summary.random$year[,2])
    params <- cbind(params, outINLA$summary.random$yearW1[,2])
    params <- cbind(params, outINLA$summary.random$yearW2[,2])
    params <- cbind(params, outINLA$summary.random$yearW3[,2])
    ngroups <- length(outINLA$summary.random$GroupID[,2])
    params <- cbind(params, c(outINLA$summary.random$GroupID[,2],rep(NA,nrow(params)-ngroups)))
    names(params) <- c("Year", "Intercept.yr", "logarea.t0.yr", 
                       "W1.yr", "W2.yr", "W3.yr", "Group")
  }
  if(ncol(crowd)==4){
    params <- as.data.frame(outINLA$summary.random$yearID[,1:2])
    params <- cbind(params, outINLA$summary.random$year[,2])
    params <- cbind(params, outINLA$summary.random$yearW1[,2])
    params <- cbind(params, outINLA$summary.random$yearW2[,2])
    params <- cbind(params, outINLA$summary.random$yearW3[,2])
    params <- cbind(params, outINLA$summary.random$yearW4[,2])
    ngroups <- length(outINLA$summary.random$GroupID[,2])
    params <- cbind(params, c(outINLA$summary.random$GroupID[,2],rep(NA,nrow(params)-ngroups)))
    names(params) <- c("Year", "Intercept.yr", "logarea.t0.yr", 
                       "W1.yr", "W2.yr","W3.yr","W4.yr", "Group")
  }
 
  #fixed effects
  fixed <- as.data.frame(outINLA$summary.fixed)[,1:2]
  tmp=matrix(NA,dim(params)[1],nrow(fixed))
  colnames(tmp)=rownames(fixed)
  tmp[1,]=fixed[,1]
  params=cbind(params,tmp)
  params$alpha=NA; params$alpha[1:length(alpha)]=alpha
  colnames(params)[which(colnames(params)=="(Intercept)")] <- "Intercept"
  
  return(list(params=params, rand_eff_summaries=rand_eff_summaries))
}






