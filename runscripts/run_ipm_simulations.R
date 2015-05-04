##  Script to run IPM simulations for each site.

rm(list=ls(all=TRUE))
library(communitySynchrony)

do_env_stoch_vec <- c(TRUE,TRUE,FALSE)
do_demo_stoch_vec <- c(TRUE,FALSE,TRUE)
sim_names <- c("ENVDEMO", "ENV", "DEMO")

Gpars_all <- readRDS("../results/growth_params_list.RDS")
Spars_all <- readRDS("../results/surv_params_list.RDS")

site_list <- names(Gpars_all)

output_list <- list()
for(do_site in site_list){
  Gpars_tmp <- Gpars_all[[do_site]]
  Spars_tmp <- Spars_all[[do_site]]
  Rpars_tmp <- Rpars_all[[do_site]]
  
  if(do_site=="Idaho"){
    Gpars_tmp <- Gpars_tmp[[-which(names(Gpars_tmp)=="ARTR")]]
    Spars_tmp <- Spars_tmp[[-which(names(Spars_tmp)=="ARTR")]]
    Rpars_tmp <- Rpars_tmp[[-which(names(Rpars_tmp)=="ARTR")]]
  }
  
  spp_list <- names(Gpars_tmp)
  Nyrs <- nrow(Gpars_tmp[[1]])
  
  # Set iteration matrix dimensions and max genet sizes by site
  # these are all taken from Chu and Adler 2015 (Ecological Monographs)
  if(do_site=="Arizona"){
    iter_matrix_dims <- c(50,50)
    max_size <- c(170,40)
  }
  if(do_site=="Idaho"){
    iter_matrix_dims <- c(75,50,50)
    max_size <- c(202,260,225)
  }
  if(do_site=="Kansas"){
    iter_matrix_dims <- c(75,50,75)
    max_size <- c(1656,823,2056)
  }
  if(do_site=="Montana"){
    iter_matrix_dims <- c(75,50,5,50)
    max_size <- c(2500,130,22,100)
  }
  if(do_site=="NewMexico"){
    iter_matrix_dims <- c(50,50)
    max_size <- c(600,1300)
  }
  
  for(i in 1:length(sim_names)){
    do_env_stoch <- do_env_stoch_vec[i]
    do_demo_stoch <- do_demo_stoch_vec[i]
    
    cover_sims <- run_ipm(A=10000, tlimit=2500, burn_in=500, spp_list=spp_list,
                          Nyrs=Nyrs, constant=do_env_stoch,
                          iter_matrix_dims=iter_matrix_dims, max_size=max_size,
                          Rpars=Rpars, Spars=Spars, Gpars=Gpars,
                          demogrpahic_stochasticity=do_demo_stoch)
    stoch_results[[sim_names[i]]] <- cover_sims
  } # end stochasticity loop
  output_list[[do_site]] <- stoch_results
} # end site loop

# Save the output
saveRDS(output_list, "../results/ipm_simulation_lists.RDS")

