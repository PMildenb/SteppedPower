% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/construct_CovMat.R
\name{construct_CovSubMat}
\alias{construct_CovSubMat}
\title{Construct a Block of the Covariance Matrix}
\usage{
construct_CovSubMat(
  N,
  timepoints,
  sigma,
  tau,
  eta = NULL,
  AR = NULL,
  rho = NULL,
  gamma = NULL,
  psi = NULL,
  INDIV_LVL = FALSE
)
}
\arguments{
\item{N}{Number of individuals per cluster}

\item{timepoints}{numeric (scalar or vector), number of timepoints (periods).
If design is swd, timepoints defaults to length(Cl)+1.
Defaults to 1 for parallel designs.}

\item{sigma}{numeric (vector of length \code{timepoints}),
residual error}

\item{tau}{numeric (vector of length \code{timepoints}),
standard deviation of random intercepts}

\item{eta}{numeric (vector of length \code{timepoints}),
standard deviation of random slope}

\item{AR}{numeric, vector containing up to three values, each between 0 and 1.
Defaults to NULL. It defines the AR(1)-correlation of random effects.
The first element corresponds to the cluster intercept, the second to the
treatment effect and the third to subject specific intercept.
If only one element is provided, autocorrelation of all random effects is
assumed to be the same.
\emph{Currently not compatible with \code{rho}!=0 !}}

\item{rho}{numeric (scalar), correlation of \code{tau} and \code{eta}. The default is no correlation.}

\item{gamma}{numeric (vector of length \code{timepoints}),
standard deviation of a random time effect.}

\item{psi}{numeric (scalar), random subject specific intercept.
Leads to a closed cohort setting}

\item{INDIV_LVL}{logical, should the computation be conducted on an
individual level? This leads to longer run time and is
mainly for diagnostic purposes.}
}
\value{
a block of a covariance matrix with two levels of clustering,
corresponding to intra-cluster covariance over time for one cluster
}
\description{
Constructs the covariance matrix
for multiple measurements of the same cluster
if the same individuals are observed at all time periods.
This function is not designed to be used directly.
}
