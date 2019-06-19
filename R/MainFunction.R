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
#' @param trt_delay numeric (possibly vector), value(s) between 0 and 1 specifing the
#' intervention effect in the first (second ... ) intervention phase
#' @param time_adjust character, specifies adjustment for time periods. Defaults to "factor".
#' @param design character, defines the type of design. Options are "SWD" and "parallel", defaults to "SWD".
#' @param EffSize numeric, raw effect
#' @param sigma numeric, residual error of cluster means if no N given.
#' @param tau numeric, standard deviation of random intercepts
#' @param eta numeric, standard deviation of random slopes **not implemented**
#' @param rho numeric, correlation of tau and eta **not implemented**
#' @param N integer, number of individuals per cluster.
#' @param Power numeric, a specified target power. If supplied, the minimal N is returned.
#' @param N_range numeric, vector specifiing the lower and upper bound for N, ignored if Power is NULL.
#' @param sig.level numeric, significance level, defaults to 0.05
#' @param df_adjust character, one of the following: **not implemented**
#' @param verbose logical, should the function return the design and covariance matrix?
#'
#' @return a list. First element is the power,
#' second element the weights of cluster per period for estimating the treatment effect **needs update**
#' @export
#'
#' @examples
#' wlsMixedPower(EffSize=1,Cl=c(1,1,1,1,1),sigma=2 ,        tau=0.2, N=c(1,1,1,1,1) )
#' wlsMixedPower(EffSize=1,Cl=c(1,1,1,1,1),sigma=2*sqrt(2) ,tau=0.2, N=c(2,2,2,2,2) )

wlsMixedPower <- function(Cl=NULL,
                          timepoints=NULL,
                          DesMat=NULL,
                          trt_delay=NULL,
                          time_adjust="factor",
                          design="SWD",
                          EffSize,sigma,
                          tau,
                          eta=NULL,
                          rho=NULL,
                          CovBlk=NULL,
                          N=NULL,
                          Power=NULL,
                          N_range=c(1,1000),
                          sig.level=0.05,
                          df_adjust="none",
                          verbose=FALSE){

  if(!is.null(N) & !is.null(Power)) stop("Both target power and individuals per cluster not NULL.")

  if(is.null(DesMat)){
    DesMat    <- construct_DesMat(Cl=Cl,trt_delay=trt_delay,design=design,
                                    timepoints=timepoints,time_adjust=time_adjust)
  }
  else if(inherits(DesMat,"list")){
    if(!inherits(DesMat[[1]],"matrix") |
       !inherits(DesMat[[2]],"numeric")|
       !inherits(DesMat[[3]],"numeric")) stop("In wlsMixedPower: Cannot interpret input for DesMat.")
    dsnmatrix     <- DesMat[[1]]
    timepoints    <- DesMat[[2]]
    SumCl         <- DesMat[[3]]
    names(DesMat) <- c("matrix","timepoints","SumCl")
  }
  else if(inherits(DesMat,"matrix")){
    dsnmatrix   <- DesMat
    timepoints  <- dim(DesMat)[2]-1
    SumCl       <- dim(DesMat)[1]/timepoints
    DesMat      <- list(matrix=dsnmatrix, timepoints=timepoints, SumCl=SumCl)
  }

  if(is.null(Power)){
    out <- wlsInnerFunction(DesMat=DesMat,
                            EffSize=EffSize,
                            sigma=sigma,
                            tau=tau,
                            Power=NULL,
                            N=N,
                            sig.level=sig.level,
                            df_adjust=df_adjust,
                            CovBlk=CovBlk,
                            verbose=verbose)
  }
  else if(Power<0 | Power>1){
    stop("Power needs to be between 0 and 1.")
  }
  else {
    N_opt <- tryCatch(ceiling(uniroot(optFunction,
                                      DesMat=DesMat,
                                      EffSize=EffSize,
                                      sigma=1,  ## why did i set this to 1 ??
                                      tau=tau,
                                      Power=Power,
                                      df_adjust=df_adjust,
                                      sig.level=.05,
                                      interval=N_range)$root),
                      error=function(cond){
                                     message(paste0("Maximal N yields power below ",Power,
                                                    ". Increase argument N_range."))
                                     return(N_range[2])})
    out <- wlsInnerFunction(DesMat=DesMat,
                            EffSize=EffSize,
                            sigma=sigma,
                            tau=tau,
                            N=N_opt,
                            sig.level=sig.level,
                            df_adjust=df_adjust,
                            CovBlk=CovBlk,
                            verbose=verbose)
    out <- append(list(N=N_opt),out)
  }
  return(out)
}

#' optFunction
#'
#' only to be called by wlsMixedPower.
#'
#' @param DesMat  list, containing a matrix, the design matrix,
#' numeric timepoints, numeric total number of Clusters
#' @param EffSize  numeric, raw effect
#' @param sigma numeric, residual error of cluster means if no N given.
#' @param tau numeric, standard deviation of random intercepts
#' @param N integer, number of individuals per cluster.
#' @param Power numeric, a specified target power. If supplied, the minimal N is returned.
#' @param sig.level numeric, significance level, defaults to 0.05

optFunction <- function(DesMat,EffSize,sigma,tau,N,Power,df_adjust,sig.level){

  diff <- (Power - wlsInnerFunction(DesMat=DesMat,
                                    EffSize=EffSize,
                                    sigma=sigma,
                                    tau=tau,
                                    N=N,
                                    df_adjust=df_adjust,
                                    sig.level=sig.level,
                                    CovBlk=NULL,
                                    verbose=F)$Power)
  return(diff)}



#' wlsInnerFunction
#'
#' @param DesMat  list, containing a matrix, the design matrix,
#' numeric timepoints, numeric total number of Clusters
#' @param EffSize  numeric, raw effect
#' @param sigma numeric, residual error of cluster means if no N given.
#' @param tau numeric, standard deviation of random intercepts
#' @param N integer, number of individuals per cluster.
#' @param Power numeric, a specified target power. If supplied, the minimal N is returned.
#' @param df_adjust character, one of the following: **not implemented**
#' @param sig.level numeric, significance level, defaults to 0.05
#' @param verbose logical, should the function return the design and covariance matrix?
#'
#' @export

wlsInnerFunction <- function(DesMat,
                             EffSize,
                             sigma,
                             tau,
                             N,
                             Power,
                             CovBlk=NULL,
                             df_adjust=df_adjust,
                             sig.level,
                             verbose){
  dsnmatrix  <- DesMat$matrix
  timepoints <- DesMat$timepoints
  SumCl      <- DesMat$SumCl
  CovMat     <- construct_CovMat(SumCl=SumCl, timepoints=timepoints,
                                    sigma=sigma, tau=tau, N=N, CovBlk=CovBlk)

  tmpmat <- t(dsnmatrix) %*% Matrix::solve(CovMat)
  VarMat <- Matrix::solve(tmpmat %*% dsnmatrix)
  WgtMat <- matrix((VarMat %*% tmpmat)[1,], nrow=SumCl, byrow=TRUE)

  df <- switch(df_adjust,
               "none"           = Inf,
               "between-within" = SumCl - rankMatrix(dsnmatrix),
               "containment"    = dim(dsnmatrix)[1] - SumCl,
               "residual"       = dim(dsnmatrix)[1] - rankMatrix(dsnmatrix))
  if(df<3){
    df <- Inf
    warning(paste0(df_adjust,"-method not applicable. No DDF adjustment used."))
  }

  out <- list(Power=tTestPwr(d=EffSize, se=sqrt(VarMat[1,1]), df=df, sig.level=sig.level))
  out <- append(out, list(denomDF=df,
                          df_adjust=df_adjust,
                          sig.level=sig.level,
                          WeightMatrix=WgtMat,
                          DesignMatrix=dsnmatrix,
                          CovarianceMatrix=CovMat))
  class(out) <- append(class(out),"wlsPower")
  return(out)
}


## Methods for swdPower

print.swdPower <- function(wlsPower){
  cat("Power                                = ", wlsPower$Power,    "\n")
  cat("ddf adjustment                       = ", wlsPower$df_adjust,"\n")
  cat("Denominator degrees of freedom       = ", wlsPower$denomDF,  "\n")
  cat("Significance level (two sided)       = ", wlsPower$sig.level,"\n")
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
#'s
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
