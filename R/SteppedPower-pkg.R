#' SteppedPower
#'
#' The package `SteppedPower` aims at providing tools for power - and sample
#' size - calculation and design diagnostics in interventional study designs,
#' with a special focus on stepped wedge designs.
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
