% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_em_data.R
\name{simulate_em_data}
\alias{simulate_em_data}
\title{Simulate OM Observation Data for TMB Assessment Model}
\usage{
simulate_em_data(
  nyears,
  dem_params,
  land_caa,
  survey_indices,
  fxfish_caa_obs,
  twfish_caa_obs,
  ll_ac_obs,
  tw_ac_obs,
  ll_srv_indic,
  tw_srv_indic,
  model_options,
  added_years = 1,
  file_suffix = ""
)
}
\arguments{
\item{nyears}{the number of years of data that will be passed to the model}

\item{dem_params}{demographic parameter matrices subsetted to 1 year}

\item{land_caa}{nyears worth of landed catch-at-age data (dimensins [1, nages, nsexes, nregions, nfleets])}

\item{survey_indices}{nyears worth of survey indices (LL RPN, LL RPW, and TW RPW)}

\item{fxfish_caa_obs}{nyears worth of catch-at-age observation from the fixed gear fishery}

\item{twfish_caa_obs}{nyears worth of catch-at-age observation from the trawl gear fishery}

\item{ll_ac_obs}{nyears worth of age composition observations frmo the longline survey}

\item{tw_ac_obs}{nyears worth of age composition observations frmo the trawl survey}

\item{model_options}{list of model options provided to the OM}

\item{file_suffix}{suffix to append to saved outputs}

\item{added_year}{the number of new years of data being added (should usually be 1)}
}
\description{
Format data and observations from an `afscOM` operating model
and update the model data and parameters required by the 
`SpatialSablefishAssessment` TMB assessment model.
}
\details{
This function requires the existence of 'data/sablefish_em_data_2022.RDS',
and 'data/sablefish_em_par_2022.RDS' files, which come packaged with the
`SablefishMSE` codebase.
}
