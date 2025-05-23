% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{aav}
\alias{aav}
\title{Average Annual Variation (AAV)}
\usage{
aav(data)
}
\arguments{
\item{data}{ordered vector of observations of a quantity}
}
\description{
Calculates average annual variation ($\frac{\sum{\frac{|x_y - x_{y-1}|}{mean(x)}}}{N-1}$)
}
