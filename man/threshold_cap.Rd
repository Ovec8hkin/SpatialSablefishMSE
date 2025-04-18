% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/harvest_control_rules.R
\name{threshold_cap}
\alias{threshold_cap}
\title{Compute F following a threshold rule, with a cap at B/B_ref = 1}
\usage{
threshold_cap(ssb, B_ref, F_ref, alpha, beta = 1)
}
\arguments{
\item{ssb}{Value for SSB}

\item{B_ref}{Value for biomass based reference point (e.g., B40)}

\item{F_ref}{Value for fishing mortality based reference point (e.g., F40)}

\item{alpha}{Cutoff to set F at 0}
}
\value{
fishing mortality rate
}
\description{
Compute F following a threshold rule, with a cap at B/B_ref = 1
}
