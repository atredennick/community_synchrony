##  Get sequence of random year effects for each site

tlimit=2000

# Read in growth parameters to get number of random years
Gpars_all <- readRDS("../results/growth_params_list.RDS")
site_list <- names(Gpars_all)

randyr_list <- list()
for(do_site in site_list){
  Gpars_tmp <- Gpars_all[[do_site]]
  Nyrs <- nrow(Gpars_tmp[[1]])
  randyr_list[[do_site]] <- sample(1:Nyrs, size = tlimit, replace = TRUE)
}

saveRDS(randyr_list, "../results/randyr_list.RDS")
