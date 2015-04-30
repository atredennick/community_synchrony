##  Functions to estimate crowding effects on focal species

#' Estimate growth crowding effects on focal species
#' 
#' @author Andrew Tredennick
#' @param site Name of the focal site (character).
#' @param data_path The directory path leading to the data folder (character).
#' @param alphas The alpha values for each species at the site, in alphabetical
#'               order by species.
#' @return List of crowding matrices, one list entry per species.

estimate_crowding <- function(site, data_path, alphas){
  species_list <- sort(list.files(data_path))
  alpha.effect <- alphas
}