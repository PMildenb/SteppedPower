
#' zTestPwr
#'
#' computes the power of a z-test given a standard error, an effect size and a significance level.
#' Computes the exact power, see second example
#'
#' @param d numeric, raw effect
#' @param se numeric, standard error
#' @param sig.level numeric, significance level, defaults to 0.05
#'
#' @return a scalar
#' @export
#'
#' @examples zTestPwr(4,1) :
#' zTestPwr(qnorm(.975),1) == pnorm(qnorm(.975),qnorm(.975),1) + pnorm(-qnorm(.975),qnorm(.975),1)


zTestPwr <- function(d,se,sig.level=0.05){
  dsz <- abs(d/se)
  Pwr <- pnorm(dsz+qnorm(sig.level/2)) + pnorm(-dsz+qnorm(sig.level/2))
  return(Pwr)
}


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
#' @examples tTestPwr(4,1,10) ; tTestPwr(4,1,30) ; zTestPwr(4,1)

tTestPwr <- function(d,se,df,sig.level=0.05){
  dsz <- abs(d/se)
  Pwr <- pt(dsz+qt(sig.level/2,df=df),df=df) + pt(-dsz+qt(sig.level/2,df=df),df=df)
  return(Pwr)
}





#' zTestSampSize
#'
#' calculate needed sample size for given target power, effect size and individual variance
#' does it work ? i fear not (got to think about it ...)
#'
#' @param d numeric, raw effect
#' @param sd numeric, standard deviaton (on individual level)
#' @param sig.level numeric, significance level, defaults to 0.05
#' @param Power target power
#'
#' @return a scalar
#' @export
#'
#' @examples zTestPwr(4,1)

zTestSampSize <- function(d,sd,Power,sig.level=0.05){
  nInd <- ((qnorm(1-sig.level/2)+qnorm(Power))*sd/d)^2
  return(nInd)
}
