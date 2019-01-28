#' wlsMixedPower
#'
#' This is the work-horse function of the SteppedPower package.
#' It calls the constructor functions for the design matrix $X$ and covariance matrix $V$,
#' and then calculates the variance of the intervention effect via
#' $$Var(trteff)=(X'V^{-1}X)^{-1}[1,1]$$
#'
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD),
#' or number in control and intervention (in parallel designs)
#' Else residual error on individual level
#' @param timepoints numeric (scalar or vector), number of timepoints (periods).
#' If design is swd, timepoints defaults to length(Cl)+1. Defaults to 1 for parallel designs.
#' @param DesMat matrix of dimension ... , if supplied, timepoints and Cl are ignored.
#' @param delay numeric (possibly vector), value(s) between 0 and 1 specifing the
#' intervention effect in the first (second ... ) intervention phase  *not implemented*
#' @param design character, defines the type of design. Options are "SWD" and "parallel", defaults to "SWD".
#' @param EffSize numeric, raw effect
#' @param sigma numeric, residual error of cluster means if no N given.
#' @param tau numeric, standard deviation of random intercepts
#' @param eta numeric, standard deviation of random slopes **not implemented**
#' @param rho numeric, correlation of tau and eta **not implemented**
#' @param N Integer,
#' @param Power numeric, a specified target power. If supplied, the minimal N is returned **work in progress**
#' @param N_range numeric, vector specifiing the lower and upper bound for N, ignored if Power is NULL.
#' @param sig.level numeric, significance level, defaults to 0.05
#' @param verbose logical, should the function return the design and covariance matrix?
#'
#' @return a list. First element is the power,
#' second element the weights of cluster per period for estimating the treatment effect **needs update**
#' @export
#'
#' @examples
#' wlsMixedPower(EffSize=1,Cl=c(1,1,1,1,1),sigma=2 ,        tau=0.2, N=c(1,1,1,1,1) )
#' wlsMixedPower(EffSize=1,Cl=c(1,1,1,1,1),sigma=2*sqrt(2) ,tau=0.2, N=c(2,2,2,2,2) )

wlsMixedPower <- function(Cl=NULL,timepoints=NULL,DesMat=NULL,delay=NULL,design="SWD",
                          EffSize,sigma,tau,eta=NULL,rho=NULL,
                          N=NULL,Power=NULL,N_range=c(1,1000),sig.level=0.05,verbose=FALSE){

  if(!is.null(N) & !is.null(Power)) stop("Both target power and individuals per cluster not NULL.")

  if(is.null(DesMat)){
    SumCl       <- sum(Cl)
    if(is.null(timepoints)){
      if(design=="SWD"){      timepoints  <- length(Cl)+1 } else
      if(design=="parallel"){ timepoints  <- 1            }
    }
    DesMat      <- construct_DesMat(Cl=Cl,delay=delay,design=design,timepoints=timepoints)
  } else {
    # try(
    #   timepoints  <- DesMat$timepoints
    #   SumCl       <- dim(DesMat)[1]/timepoints
    # ) else {
      timepoints  <- dim(DesMat[2])-1                 ## this is only the 'canonical' case
      SumCl       <- dim(DesMat)[1]/timepoints
    # }

  }

  if(is.null(Power)){
    out <- wlsInnerFunction(DesMat=DesMat,EffSize=EffSize,SumCl=SumCl,timepoints=timepoints,
                            sigma=sigma,tau=tau,Power=NULL,N=N,sig.level=sig.level,
                            verbose=verbose)
  }else{
    optFunction <- function(DesMat,EffSize,SumCl,timepoints,sigma,tau,Power,N,sig.level,verbose){

      diff <- abs(Power - wlsInnerFunction(DesMat=DesMat,SumCl=SumCl,timepoints=timepoints,
                                           EffSize=EffSize,sigma=sigma,tau=tau,
                                           N=N,sig.level=sig.level,verbose=verbose)$`Power`)
      return(diff)}

    N_opt <- ceiling(optim(par=sqrt(N_range[2]),optFunction,DesMat=DesMat,EffSize=EffSize,SumCl=SumCl,
                           timepoints=timepoints,sigma=1,tau=tau,
                           Power=Power,sig.level=.05,verbose=F,
                           method="Brent",lower=N_range[1],upper=N_range[2])$`par`)

    out <- wlsInnerFunction(DesMat=DesMat,EffSize=EffSize,SumCl=SumCl,timepoints=timepoints,
                     sigma=sigma,tau=tau,N=N_opt,sig.level=sig.level,
                     verbose=verbose)
    out <- append(list(N=N_opt),out)
  }
  return(out)
}


#' wlsInnerFunction
#'
#' only to be called by wlsMixedPower
#'
#' @param DesMat matrix, the design matrix
#' @param SumCl integer, total number of clusters
#'
#' @return a list.
#'
#' @export
#'

wlsInnerFunction <- function(DesMat,EffSize,SumCl,timepoints,sigma,tau,N,
                             Power,sig.level,verbose){

  CovMat        <- construct_CovMat(SumCl=SumCl,timepoints=timepoints,
                                    sigma=sigma,tau=tau,N=N)

  tmpmat <- t(DesMat) %*% Matrix::solve(CovMat)
  VarMat <- Matrix::solve(tmpmat %*% DesMat)
  WgtMat <- matrix((VarMat %*% tmpmat)[1,],nrow = SumCl,byrow=TRUE)

  out <- list(Power=zTestPwr(d=EffSize,se=sqrt(VarMat[1,1]),sig.level=sig.level))
  if(verbose)
    out <- append(out, list(WeightMatrix=WgtMat, DesignMatrix=DesMat, CovarianceMatrix=CovMat))
  return(out)
}



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



#' zTestSampSize
#'s
#' calculate needed sample size for given target power, effect size and individual variance
#' does it work ? i fear not (gotta think about it ...)
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
