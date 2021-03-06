% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/run_ipm.R
\name{run_ipm}
\alias{run_ipm}
\title{Run multi-species Integral Projection Model}
\usage{
run_ipm(A = 10000, tlimit = 2500, burn_in = 500, spp_list, Nyrs,
  constant = FALSE, iter_matrix_dims, max_size, Rpars, Spars, Gpars,
  demographic_stochasticity = FALSE)
}
\arguments{
\item{A}{Area of quadrat (cm X cm).}

\item{tlimit}{Length of simulation.}

\item{burn_in}{Number of years to discard before calculating things.}

\item{spp_list}{Character vector with list of species for focal site.}

\item{Nyrs}{Number of years represented in the yearly regression coefficients.}

\item{constant}{TRUE/FALSE flag for simulating in a constant environment (no random year effects).}

\item{iter_matrix_dims}{The size of the big matrix for each species.}

\item{max_size}{Maximum size, in cm^2, for allowable for each species (in alphabetical order).}

\item{Rpars}{Recruitment regression parameters for focal site and species.}

\item{Spars}{Survival regression parameters for focal site and species.}

\item{Gpars}{Recruitment regression parameters for focal site and species.}

\item{demographic_stochasticity}{TRUE/FALSE flag for including demographic stoachasticity or not.}
}
\description{
Run multi-species Integral Projection Model
}
\author{
Andrew Tredennick
}

