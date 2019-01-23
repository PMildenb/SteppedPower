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
#' @param verbose logical, should the function return the design and covariance matrix?
#'
#' @return a list. First element is the power,
#' second element the weights of cluster per period for estimating the treatment effect
#' @export
#'
#' @examples
#' wlsMixedPower(EffSize=1,Cl=c(1,1,1,1,1),sigma=2 ,        tau=0.2, N=c(1,1,1,1,1) )
#' wlsMixedPower(EffSize=1,Cl=c(1,1,1,1,1),sigma=2*sqrt(2) ,tau=0.2, N=c(2,2,2,2,2) )

wlsMixedPower <- function(Cl=NULL,timepoints=NULL,delay=NULL,
                          DesMat=NULL,design="SWD",family="gaussian",
                          EffSize,sigma,tau,
                          N=NULL,Power=NULL,sig.level=0.05,verbose=FALSE){

  if(!is.null(N) & !is.null(Power)) stop("Both target power and individuals per cluster not NULL.")

  if(is.null(DesMat)){
    SumCl       <- sum(Cl)
    if(is.null(timepoints)){
      if(design=="SWD"){      timepoints  <- length(Cl)+1 } else
      if(design=="parallel"){ timepoints  <- 1            }
    }
    DesMat      <- construct_DesMat(Cl=Cl,delay=delay,design=design,timepoints=timepoints)
  } else {
    timepoints  <- dim(DesMat)[2]-1                 ## this is only the 'canonical' case
    SumCl       <- dim(DesMat)[1]/timepoints
  }

  if(is.null(Power)){
    out <- wlsInnerFunction(DesMat=DesMat,EffSize=EffSize,SumCl=SumCl,timepoints=timepoints,
                            sigma=sigma,tau=tau,family=family,Power=NULL,N=N,sig.level=sig.level,
                            verbose=verbose)
  }else{
    optFunction <- function(DesMat,EffSize,SumCl,timepoints,sigma,tau,family,Power,N,sig.level,verbose){

      diff <- abs(Power - wlsInnerFunction(DesMat=DesMat,EffSize=EffSize,SumCl=SumCl,timepoints=timepoints,
                               sigma=sigma,tau=tau,family=family,N=N,sig.level=sig.level,
                               verbose=verbose)$`Power`)
      return(diff)}

    N_opt <- ceiling(optim(par=1,optFunction,DesMat=DesMat,EffSize=EffSize,SumCl=SumCl,
                           timepoints=timepoints,sigma=1,tau=.3,family="gaussian",
                           Power=.9,sig.level=.05,verbose=F,
                           method="Brent",lower=1,upper=1000)$`par`)

    out <- wlsInnerFunction(DesMat=DesMat,EffSize=EffSize,SumCl=SumCl,timepoints=timepoints,
                     sigma=sigma,tau=tau,family=family,N=N_opt,sig.level=sig.level,
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

wlsInnerFunction <- function(DesMat,EffSize,SumCl,timepoints,sigma,tau,family,N,
                             Power,sig.level,verbose){

  CovMat        <- construct_CovMat(SumCl=SumCl,timepoints=timepoints,
                                    sigma=sigma,tau=tau,family=family,
                                    N=N)

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
#' @param se numeric, standard deviaton (on individual level)
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
