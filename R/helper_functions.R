
#' tTestPwr
#'
#' computes the power of a t-test given a standard error, an effect size,
#' the degrees of freedom of the t-distribution and a significance level.
#' Computes the exact power, see second example
#'
#' @param d numeric, raw effect
#' @param se numeric, standard error
#' @param df numeric, degrees of freedom of the t-distribution
#' @param sig.level numeric, significance level, defaults to 0.05
#'
#' @return a scalar
#' @export
#'
#' @examples tTestPwr(4,1,10) ; tTestPwr(4,1,30) ; tTestPwr(4,1,Inf)

tTestPwr <- function(d,se,df,sig.level=0.05){
  dsz <- abs(d/se)
  q   <- qt(sig.level/2,df=df)
  Pwr <- pt(dsz + q, df=df) + pt(-dsz + q, df=df)
  return(Pwr)
}

# tTestPwr uses a central t-distribution. As applied sometimes in the literature,
# e.g. Li, Turner, Preisser 2018 | Li, Redden 2015
# This is NOT the same as a t-Test, which uses a non-central t-distribution !!
# Hence the name of the above function should be changed to 'scaledWaldPwr'
# tTestPwr2 tries to mimic a t-test when only df are given.

tTestPwr2 <- function(d,se,df,sig.level=.05){
  d   <- abs(d/se)
  q   <- qt(sig.level/2, df=df)
  ncp <- -sqrt((df+2)/4)*d ## TODO: ?? how to calculate n using df ??
  pwr <- pt(q,  df, ncp=ncp) + pt(-q, df,ncp=ncp,lower=FALSE)
  return(pwr)
}
