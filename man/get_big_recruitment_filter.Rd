% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_big_recruitment_filter}
\alias{get_big_recruitment_filter}
\title{Create dplyr Filter for Years Around Large Recruitment Events}
\usage{
get_big_recruitment_filter(
  om_list,
  seed_list,
  model_runs,
  om_names,
  large_event_thresh,
  lags
)
}
\arguments{
\item{om_list}{a list of OM objects (as would be given to `run_mse_multiple`)}

\item{seed_list}{a vector or list of simulation seeds (as would be given to `run_mse_multiple`)}

\item{model_runs}{list of completed MSE model runs (the output of `run_mse_multiple`)}

\item{large_event_thresh}{threshold value for declaring a recruitment event "large"}
}
\description{
Given a list of MSE model runs, identifies years in which recruitment was
larger than a specified threshold value for each OM simulation. That information
is used to generate an complex boolean expression object for use with 
dplyr::filter(...) when computing performance metrics.
}
