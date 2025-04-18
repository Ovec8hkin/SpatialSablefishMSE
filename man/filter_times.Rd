% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_utils.R
\name{filter_times}
\alias{filter_times}
\title{Filter processed MSE data to between two time points,}
\usage{
filter_times(data, time_horizon)
}
\arguments{
\item{data}{dataframe of processed MSE data}

\item{time_horizon}{vector of times to filter between.
If first element is NA, lower bound will be 1. If
second element is NA, upper bound will be maximum value
in `time` column of `data`.}
}
\description{
Filter dataframe of processed MSE data by 'time'
column to between times set by `time_horizon`.
}
