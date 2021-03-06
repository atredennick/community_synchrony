% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/recruitment_function.R
\name{recruit_mcmc}
\alias{recruit_mcmc}
\title{Estimate recruitment regression coefficients using MCMC}
\usage{
recruit_mcmc(dataframe, n_adapt = 5000, n_update = 10000,
  n_samples = 20000, n_thin = 50, sppList)
}
\arguments{
\item{dataframe}{Merged time series dataframe of recruitment area for all species.}

\item{n_adapt}{Number of iterations for adaptation of MCMC (default = 5000).}

\item{n_update}{Number of iterations for update phase of MCMC (default = 10000).}

\item{n_samples}{Number of samples to collect from MCMC after update phase (default = 20000).}

\item{n_thin}{Number of iterations by which to thin the samples (default = 50).}

\item{sppList}{Character vector of species code names for the site.}
}
\value{
Matrix of statistical results by fitted parameter.
}
\description{
Estimate recruitment regression coefficients using MCMC
}
\author{
Andrew Tredennick
}

