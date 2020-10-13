#' SteppedPower
#'
#' SteppedPower offers options for power calculation as well
#' as sample size calculation (by the use of uniroot() ) for
#' longitudinal mixed model settings with up to two levels of clustering.
#' All calculations are oracle estimates i.e. assume random effect variances
#' to be known (or guessed) in advance.
#' @import knitr
#' @importFrom stats coef family gaussian optim pnorm qnorm dnorm rbinom rnorm
#' uniroot integrate pt qt binomial
#' @import Matrix
#' @importFrom grDevices colorRamp
#' @importFrom utils adist
#' @import plotly
#'
#'
#'
#' @author Philipp Mildenberger \email{pmildenb@@uni-mainz.de}
#' @name SteppedPower-pkg
#' @docType package
NULL
