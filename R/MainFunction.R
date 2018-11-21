#'
#' 'wlsMixedPower'
#'
#' @param EffSize numeric, raw effect
#' @param I integer (vector), number of clusters per wave (in SWD)
#' @param sigma numeric, residual error of cluster means if no N given.
#' Else residual error on individual level
#' @param tau numeric, standard deviation of random intercepts
#' @param sig.level numeric, significance level, defaults to 0.05
#' @param delay numeric (possibly vector), value between 0 and 1 specifing the
#' intervention effect in the first (second ... ) intervention phase  *not implemented*
#'
#' @return a formula for CoxEval()
#' @export
#'
#' @examples
#' BuildFormula(VarName="Variable",AdjName="additiveAdjustment",StratName="Stratification")

wlsMixedPower <- function(EffSize,sigma,tau,family=gaussian(),timepoints=NULL,
                          N=NULL,sig.level=0.05,DesMat=NULL,I=NULL,delay=NULL){

  if(is.null(DesMat)){
    DesMat      <- construct_DesMat(I=I,delay=delay)
    timepoints  <- length(I)+1
    SumCl       <- sum(I)
  } else {
    timepoints  <- dim(DesMat)[2]-1                 ## special case
    SumCl       <- dim(DesMat)[1]/timepoints
  }

  CovMat        <- construct_CovMat(SumCl=SumCl,timepoints=timepoints,
                                    sigma=sigma,tau=tau,family=family,
                                    N=N)

  tmpmat <- t(DesMat) %*% solve(CovMat)
  VarMat <- solve(tmpmat %*% DesMat)
  HatMat <- matrix((VarMat %*% tmpmat)[1,],nrow = SumCl,byrow=TRUE)  ## HatMatrix solle funktionieren
  Pwr    <- zTestPwr(d=EffSize,se=sqrt(VarMat[1,1]),sig.level=sig.level)
  return(list(Power=Pwr,HatMatrix=HatMat))
}


#'
#' 'zTestPower'
#'
#' @param d numeric, raw effect
#' @param se numeric, standard error
#' @param sig.level numeric, significance level, defaults to 0.05
#'
#' @return a formula for CoxEval()
#' @export
#'
#' @examples
#' BuildFormula(VarName="Variable",AdjName="additiveAdjustment",StratName="Stratification")


zTestPwr <- function(d,se,sig.level=0.05){
  dsz <- abs(d/se)
  Pwr <- pnorm(dsz+qnorm(sig.level/2)) + pnorm(-dsz+qnorm(sig.level/2))
  return(Pwr)
}
