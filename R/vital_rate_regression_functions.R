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

get_growth_params <- function(dataframe, crowd_mat, alpha){
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
  
  formula <- logarea.t1 ~ logarea.t0+crowd+
    f(yearID, model="iid", prior="normal",param=c(0,0.001))+
    f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
    f(year, logarea.t0, model="iid", prior="normal",param=c(0,0.001))
  
  outINLA <- inla(formula, data=D,
                  family=c("gaussian"), verbose=FALSE,
                  control.predictor = list(link = 1),
                  control.compute=list(dic=T,mlik=T),
                  control.inla = list(h = 1e-6))
  
  # Fit variance
  x <- outINLA$summary.fitted.values$mean #fitted values from INLA
  y <- (D$logarea.t1-outINLA$summary.fitted.values$mean)^2 #calculates the variance (residuals^2)
  outVar <- nls(y~a*exp(b*x),start=list(a=1,b=0)) #model the variance as function of fitted size
  
  # Collect parameters
  # random year and group effects
  params <- as.data.frame(outINLA$summary.random$yearID[,1:2])
  params <- cbind(params, outINLA$summary.random$year[,2])
  names(params) <- c("Year", "Intercept.yr", "logarea.t0.yr")
  tmp <- as.data.frame(outINLA$summary.random$Group[,2])
  names(tmp) <- "Group"
  tmp$GrpName <- outINLA$summary.random$Group[,1]
  tmp[(NROW(tmp)+1):NROW(params),]=NA
  params=cbind(params,tmp)
  
  # fixed effects
  fixed <- as.data.frame(outINLA$summary.fixed)[,1:2]
  tmp=matrix(NA,dim(params)[1],nrow(fixed))
  colnames(tmp)=rownames(fixed)
  tmp[1,]=fixed[,1]
  params=cbind(params,tmp)
  colnames(params)[6] <- "Intercept"
  params$alpha=NA; params$alpha[1:length(alpha)]=alpha
  #variance 
  params$sigma.a=NA; params$sigma.a[1]=coef(outVar)[1] 
  params$sigma.b=NA; params$sigma.b[1]=coef(outVar)[2]
  
  return(params)
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

get_survival_params <- function(dataframe, crowd_mat, alpha){
  D <- dataframe
  D$logarea <- log(D$area)
  D$quad <- as.character(D$quad)
  
  crowd = crowd_mat 
  crowd[crowd<1e-99]=0 # make really small crowding indices 0
  
  library(INLA)
  # Set up ID variables for INLA random effects
  D$yearID <- D$year+max(D$year) # for random year offset on intercept
  D$GroupID <- as.numeric(D$Group)
  
  formula2 <- survives ~ logarea+crowd+
    f(yearID, model="iid", prior="normal",param=c(1,0.001))+
    f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
    f(year, logarea, model="iid", prior="normal",param=c(0,0.001))
  
  outINLA <- inla(formula2, data=D,
                  family=c("binomial"), verbose=FALSE,
                  control.compute=list(dic=T,mlik=T),
                  control.predictor = list(link = 1),
                  control.inla = list(h = 1e-6),
                  Ntrials=rep(1,nrow(D)))
  
  #Collect parameters
  #random year and group effects
  params <- as.data.frame(outINLA$summary.random$yearID[,1:2])
  params <- cbind(params, outINLA$summary.random$year[,2])
  ngroups <- length(outINLA$summary.random$GroupID[,2])
  params <- cbind(params, c(outINLA$summary.random$GroupID[,2],rep(NA,nrow(params)-ngroups)))
  names(params) <- c("Year", "Intercept.yr", "logarea.yr", "Group")
 
  #fixed effects
  fixed <- as.data.frame(outINLA$summary.fixed)[,1:2]
  tmp=matrix(NA,dim(params)[1],nrow(fixed))
  colnames(tmp)=rownames(fixed)
  tmp[1,]=fixed[,1]
  params=cbind(params,tmp)
  params$alpha=NA; params$alpha[1:length(alpha)]=alpha
  colnames(params)[5] <- "Intercept"
  
  return(params)
}






