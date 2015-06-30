##  Individual-based model to test poisson approximation 
##  to demographic stochasticity.

##  We ignore crowding so it can be spatially implicit

library(communitySynchrony)

####
####  Import vital rate regression parameter functions
####

do_site <- "Idaho"

Gpars_all <- readRDS("../../results/growth_params_list.RDS")
Spars_all <- readRDS("../../results/surv_params_list.RDS")
Rpars_all <- readRDS("../../results/recruit_parameters.RDS")
Gpars_tmp <- Gpars_all[[do_site]]
Spars_tmp <- Spars_all[[do_site]]
Rpars_tmp <- Rpars_all[[do_site]]

spp_list <- names(Gpars_tmp)
Nyrs <- nrow(Gpars_tmp[[1]])

site_path <- paste("../../data/", do_site, sep="")
Gpars <- format_growth_params(do_site = do_site, species_list = spp_list, 
                              Nyrs = Nyrs, Gdata_species = Gpars_tmp)
Spars <- format_survival_params(do_site = do_site, species_list = spp_list, 
                                Nyrs = Nyrs, Sdata_species = Spars_tmp)
Rpars <- format_recruitment_params(do_site = do_site, species_list = spp_list, 
                                   Nyrs = Nyrs, Rdata_species = Rpars_tmp,
                                   path_to_site_data = site_path) 

####
####  Make vital rate functions
####
survival <- function(Spars, plants){
  
}



