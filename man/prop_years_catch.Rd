% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_metrics.R
\name{prop_years_catch}
\alias{prop_years_catch}
\title{Compute proportion of years catch exceeds a threshold level across projection period}
\usage{
prop_years_catch(
  model_runs,
  extra_columns,
  hcr_filter,
  om_filter,
  catch_threshold,
  interval_widths = c(0.5, 0.8),
  time_horizon = c(65, NA),
  extra_filter = NULL,
  relative = NULL,
  summarise_by = c("om", "hcr"),
  summary_out = TRUE
)
}
\arguments{
\item{model_runs}{list of completed MSE simulations runs}

\item{extra_columns}{data.frame specifying names for OM and HCR to attach
to each model_run (see `bind_mse_outputs` for more details)}

\item{hcr_filter}{vector of HCR names to calculate metric over}

\item{om_filter}{vector of OM names to calculate metric over}

\item{interval_widths}{confidence intevrals to compute}

\item{extra_filter}{an additional set of filters to apply before computing 
medians and confidence intervals}

\item{relative}{a management procedure to compute metric relative to}

\item{summarise_by}{vector of columns to summarise metric by}

\item{summary_out}{whether to output data summarised by `ggdist` or full data}
}
\description{
Compute average proportion of years where catch exceeds a threshold level
(median and CIs) per year across all years and simulation seeds, 
for each combination of operating models and management procedures.
}
