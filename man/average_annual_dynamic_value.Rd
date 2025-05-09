% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_metrics.R
\name{average_annual_dynamic_value}
\alias{average_annual_dynamic_value}
\title{Compute Average Annual Dynamic Value of Catch across projection period}
\usage{
average_annual_dynamic_value(
  model_runs,
  extra_columns,
  interval_widths = c(0.5, 0.8),
  time_horizon = c(65, NA),
  extra_filter = NULL,
  relative = NULL,
  summarise_by = c("om", "hcr")
)
}
\arguments{
\item{model_runs}{list of completed MSE simulations runs}

\item{extra_columns}{data.frame specifying names for OM and HCR to attach
to each model_run (see `bind_mse_outputs` for more details)}

\item{interval_widths}{confidence intevrals to compute}

\item{extra_filter}{an additional set of filters to apply before computing 
medians and confidence intervals}

\item{relative}{a management procedure to compute metric relative to}

\item{summarise_by}{vector of columns to summarise metric by}
}
\description{
Compute average annual dynamic value of landed catch from the
fixed-gear fleet (median and CIs) per year across all years and 
simulation seeds, for each combination of operating models and management 
procedures. Value is computed based on an assumed linear relationship
between catch and price (assuming that increased catch floods existing
markets and decreases prices due to market abundance), where prices are 
constant at a "maximum" level when catch <15k mt, and are constant at a 
"minimum" level when catch >30k mt. Pricing data from June 2024 NPFMC
Small Sablefish Release Analyses.
}
