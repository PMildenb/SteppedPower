
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
