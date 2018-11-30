#' wlsMixedPower
#'
#' This is the work-horse function of the SteppedPower package.
#' It calls the constructor functions for the design matrix $X$ and covariance matrix $V$,
#' and then calculates the variance of the intervention effect via
#' $$Var(trteff)=(X'V^{-1}X)^{-1}[1,1]$$
#'
#'
#' @param EffSize numeric, raw effect
#' @param Cl integer (vector), number of clusters per wave (in SWD),
#' or number in control and intervention (in parallel designs)
#' @param sigma numeric, residual error of cluster means if no N given.
#' Else residual error on individual level
#' @param tau numeric, standard deviation of random intercepts
#' @param sig.level numeric, significance level, defaults to 0.05
#' @param delay numeric (possibly vector), value between 0 and 1 specifing the
#' intervention effect in the first (second ... ) intervention phase  *not implemented*
#' @param design character, defines the type of design. Options are "SWD" and "parallel", defaults to "SWD".
#'
#' @return a list of lenght two. First element is the power,
#' second element the weights of cluster per period for estimating the treatment effect
#' @export
#'
#' @examples
#' wlsMixedPower(EffSize=1,Cl=c(1,1,1,1,1),sigma=2 ,        tau=0.2, N=c(1,1,1,1,1) )
#' wlsMixedPower(EffSize=1,Cl=c(1,1,1,1,1),sigma=2*sqrt(2) ,tau=0.2, N=c(2,2,2,2,2) )

wlsMixedPower <- function(EffSize,sigma,tau,family=gaussian(),timepoints=NULL,
                          N=NULL,sig.level=0.05,DesMat=NULL,Cl=NULL,delay=NULL,design="SWD"){
  if(is.null(DesMat)){
    SumCl       <- sum(Cl)
    if(is.null(timepoints)){
      if(design=="SWD"){      timepoints  <- length(Cl)+1 } else
      if(design=="parallel"){ timepoints  <- 1            }
    }
    DesMat      <- construct_DesMat(Cl=Cl,delay=delay,design=design,timepoints=timepoints)
  } else {
    timepoints  <- dim(DesMat)[2]-1                 ## this is only a special case
    SumCl       <- dim(DesMat)[1]/timepoints
  }

  CovMat        <- construct_CovMat(SumCl=SumCl,timepoints=timepoints,
                                    sigma=sigma,tau=tau,family=family,
                                    N=N)

  tmpmat <- t(DesMat) %*% Matrix::solve(CovMat)
  VarMat <- Matrix::solve(tmpmat %*% DesMat)
  WgtMat <- matrix((VarMat %*% tmpmat)[1,],nrow = SumCl,byrow=TRUE)
  Pwr    <- zTestPwr(d=EffSize,se=sqrt(VarMat[1,1]),sig.level=sig.level)

  return(list(Power=Pwr, WeightMatrix=WgtMat, DesignMatrix=DesMat, CovarianceMatrix=CovMat))
}


#' zTestPower
#'
#' @param d numeric, raw effect
#' @param se numeric, standard error
#' @param sig.level numeric, significance level, defaults to 0.05
#'
#' @return a scalar
#' @export
#'
#' @examples zTestPwr(4,1)


zTestPwr <- function(d,se,sig.level=0.05){
  dsz <- abs(d/se)
  Pwr <- pnorm(dsz+qnorm(sig.level/2)) + pnorm(-dsz+qnorm(sig.level/2))
  return(Pwr)
}
