#'
#' 'zTestPower'
#'
#' @param d numeric, raw effect
#' @param se numeric, standard error
#' @param sig.level numeric, significance level, defaults to 0.05
#'
#'
#'
#'

zTestPwr <- function(d,se,sig.level=0.05){
  dsz <- abs(d/se)
  Pwr <- pnorm(dsz+qnorm(sig.level/2)) + pnorm(-dsz+qnorm(sig.level/2))
  return(Pwr)
}



#'
#' 'wlsMixedPower'
#'
#' @param EffSize numeric, raw effect
#' @param I integer (vector), number of clusters per wave (in SWD)
#' @param sigma numeric, residual error of cluster means if no N given.
#' Else residual error on individual level
#' @param tau numeric, standard deviation of random intercepts
#' @param sig.level numeric, significance level, defaults to 0.05
#'
#'
#'
#'

wlsMixedPower <- function(EffSize,I,sigma,tau,family=gaussian(),N=NULL,sig.level=0.05){

  DesMat <- construct_DesMat(I=I)
  CorMat <- construct_CorMat(I=I,sigma=sigma,tau=tau,family=family,N=N)

  VarTrt <- solve(t(DesMat) %*% solve(CorMat) %*% DesMat)[1,1]

  Pwr    <- zTestPwr(d=EffSize,se=sqrt(VarTrt),sig.level=sig.level)
  return(Pwr)
}

