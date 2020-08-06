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
#' @param trt_delay numeric (possibly vector), value(s) between 0 and 1 specifying the
#' intervention effect in the first (second ... ) intervention phase
#' @param time_adjust character, specifies adjustment for time periods. Defaults to "factor".
#' @param design character, defines the type of design. Options are "SWD" and "parallel", defaults to "SWD".
#' @param EffSize numeric (scalar), raw effect
#' @param sigma numeric, residual error of cluster means if no N given.
#' @param tau numeric, standard deviation of random intercepts
#' @param eta numeric, standard deviation of random slopes
#' @param rho numeric, correlation of tau and eta **not implemented**
#' @param N numeric, number of individuals per cluster. Either a scalar, vector
#' of length #Clusters or a matrix of dimension #Clusters x timepoints
#' @param family character, distribution family. Defaults to "gaussian" **not implemented**
#' @param Power numeric, a specified target power. If supplied, the minimal N is returned.
#' @param N_range numeric, vector specifying the lower and upper bound for N, ignored if Power is NULL.
#' @param sig.level numeric, significance level, defaults to 0.05
#' @param df_adjust character, one of the following: **not implemented**
#' @param verbose logical, should the function return the design and covariance matrix?
#' @param period numeric (scalar)
#' @param CovMat numeric, a positive-semidefinite matrix of dimension (#Clusters \cdot timepoints) **experimental**
#' *#Cluster* $\cdot$ *timepoints* rows/columns. If `CovMat` is given, `sigma`,
#' `tau`, `eta` and `rho` are ignored.
#'
#' @return a list. First element is the power,
#' second element the weights of cluster per period for estimating the treatment effect **needs update**
#' @export
#'
#' @examples
#' wlsMixedPower(EffSize=1,Cl=c(1,1,1,1,1),sigma=2 ,        tau=0.2, N=c(1,1,1,1,1) )
#' wlsMixedPower(EffSize=1,Cl=c(1,1,1,1,1),sigma=2*sqrt(2) ,tau=0.2, N=c(2,2,2,2,2) )

wlsMixedPower <- function(Cl            =NULL,
                          timepoints    =NULL,
                          DesMat        =NULL,
                          trt_delay     =NULL,
                          time_adjust   ="factor",
                          period        =NULL,
                          design        ="SWD",
                          EffSize,
                          sigma         =1, ## default needed for CovMat input, still experimental
                          tau           =0,
                          eta           =NULL,
                          rho           =NULL,
                          gamma         =NULL,
                          CovMat        =NULL,
                          N             =NULL,
                          Power         =NULL,
                          family        ="gaussian",
                          N_range       =c(1,1000),
                          sig.level     =0.05,
                          df_adjust     ="none",
                          verbose       =FALSE){
  ## CHECKS #####
  if(!is.null(N) & !is.null(Power))
    stop("Both target power and individuals per cluster not NULL.")

  ## check DesMat #####
  if(is.null(DesMat)){
    DesMat    <- construct_DesMat(Cl=Cl,trt_delay=trt_delay,design=design,
                                  timepoints=timepoints,time_adjust=time_adjust,
                                  period=period)
  }else if(inherits(DesMat,"matrix") & !inherits(DesMat,"DesMat")){
    DesMat <- construct_DesMat(trtmatrix=DesMat) ## TODO : dimension checks
  }else if(!inherits(DesMat,"DesMat"))
    stop("In wlsMixedPower: Cannot interpret input for DesMat. ",
         "It must be either an object of class DsnMat or a matrix")

  ## calculate samplesize (if needed) #####
  if(!is.null(Power)){
    if(Power<0 | Power>1) stop("Power needs to be between 0 and 1.")
    N_opt <- tryCatch(ceiling(
              uniroot(function(N){Power - compute_wlsPower(DesMat    =DesMat,
                                                           EffSize   =EffSize,
                                                           sigma     =sigma,
                                                           tau       =tau,
                                                           eta       =eta,
                                                           rho       =rho,
                                                           gamma     =gamma,
                                                           N         =N,
                                                           df_adjust =df_adjust,
                                                           sig.level =sig.level,
                                                           CovMat    =CovMat,
                                                           verbose   =FALSE)$Power},
                interval=N_range)$root),
              error=function(cond){
                message(paste0("Maximal N yields power below ",Power,". Increase argument N_range."))
                return(N_range[2])
              })
    N <- N_opt
  }
  ## calculate Power #####
  out <- compute_wlsPower(DesMat    =DesMat,
                          EffSize   =EffSize,
                          sigma     =sigma,
                          tau       =tau,
                          eta       =eta,
                          rho       =rho,
                          gamma     =gamma,
                          N         =N,
                          df_adjust =df_adjust,
                          sig.level =sig.level,
                          CovMat    =CovMat,
                          verbose   =verbose)
  if(!is.null(Power)) out$N_opt <- N_opt

  return(out)
}


#' compute_wlsPower
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
#' @param rho numeric, correlation of tau and eta **not implemented**
#' @param CovMat numeric, a positive-semidefinite matrix  with
#' *#Cluster* $\cdot$ *timepoints* rows/columns. If `CovMat` is given, `sigma`,
#' `tau`, `eta` and `rho` are ignored.
#'
#' @export

compute_wlsPower <- function(DesMat,
                             EffSize,
                             sigma,
                             tau        =0,
                             eta        =NULL,
                             rho        =NULL,
                             gamma      =NULL,
                             N          =NULL,
                             CovMat     =NULL,
                             df_adjust  ="none",
                             sig.level  =.05,
                             verbose    =FALSE){
  dsnmatrix  <- DesMat$dsnmatrix
  timepoints <- DesMat$timepoints
  SumCl      <- DesMat$SumCl
  trtMat     <- DesMat$trtMat

  ## get covariance matrix #####
  if(is.null(CovMat))
    CovMat   <- construct_CovMat(SumCl      =SumCl,
                                 timepoints =timepoints,
                                 sigma      =sigma,
                                 tau        =tau,
                                 eta        =eta,
                                 rho        =rho,
                                 gamma      =gamma,
                                 trtMat     =trtMat,
                                 N          =N)

  ## matrices for power calculation #####
  tmpmat <- t(dsnmatrix) %*% Matrix::chol2inv(Matrix::chol(CovMat))
  VarMat <- Matrix::solve(tmpmat %*% dsnmatrix)
  if(verbose) ProjMat <- matrix((VarMat %*% tmpmat)[1,], nrow=SumCl, byrow=TRUE)

  ## ddf for power calculation #####
  df <- switch(df_adjust,
               "none"           = Inf,
               "between-within" = SumCl - rankMatrix(dsnmatrix),
               "containment"    = dim(dsnmatrix)[1] - SumCl,
               "residual"       = dim(dsnmatrix)[1] - rankMatrix(dsnmatrix))
  if(df<3){
    warning(df_adjust,"-method not applicable. No DDF adjustment used.")
    df <- Inf }

  Pwr <- tTestPwr(d=EffSize, se=sqrt(VarMat[1,1]), df=df, sig.level=sig.level)
  out <- list(Power     =Pwr,
              Params =list(N         =N,
                           sigma     =sigma, ## NOT compatible with CovMat-Input (!)
                           tau       =tau,
                           eta       =eta,
                           rho       =rho,
                           gamma     =gamma,
                           denomDF   =df,
                           df_adjust =df_adjust,
                           sig.level =sig.level))
  if(verbose)
    out <- append(out,
                  list(ProjMatrix     =ProjMat,
                       DesignMatrix     =DesMat,
                       CovarianceMatrix =CovMat))
  class(out) <- append(class(out),"wlsPower")
  return(out)
}

#' print.wlsPower
#'
#' @param x w
#' @param ... Arguments to be passed to methods
#'
#' @method print wlsPower
#'
#' @export
#'
#'
print.wlsPower <- function(x, ...){
  cat("Power                                = ", x$Power,    "\n")
  cat("ddf adjustment                       = ", x$Params$df_adjust,"\n")
  if(x$Params$df_adjust!="none")
  cat("Denominator degrees of freedom       = ", x$Params$denomDF,  "\n")

  cat("Significance level (two sided)       = ", x$Params$sig.level,"\n")

  if("N_opt" %in% names(x))
  cat("Needed N per cluster per period      = ", x$N_opt,"\n" )
}



#' plot.wlsPower
#'
#' @param x object of class wlsPower
#'
#' @method plot wlsPower
#'
#' @export
#'
plot.wlsPower <- function(x){
  if(!"ProjMatrix" %in% names(x)) stop("Please rerun wlsMixedPower with `verbose=TRUE` ")
  wgt <- x$ProjMatrix
  mx <- max(abs(wgt))
  sumCl <- dim(wgt)[1]
  timep <- dim(wgt)[2]
  subp <- subplot(
    plot_ly(data=data.frame(time=1:dim(wgt)[2], weight=colSums(abs(wgt))),  ## arithmetic or harmonic mean/sum ??
            type="bar", x=~time, y=~weight, color=I("grey")) %>%
      layout(xaxis=list(showticklabels=FALSE)),
    plotly_empty(type="scatter",mode="marker"),
    plot_ly(x=1:timep,y=1:sumCl,z=wgt,type="heatmap",
            colors=colorRamp(c("steelblue","white","firebrick")),
            xgap=.3,ygap=.3) %>%
      colorbar(len=1,limits=c(-mx,mx)) %>% layout(yaxis = list(autorange = "reversed")),
    plot_ly(data=data.frame(cluster=1:dim(wgt)[1], weight=rowSums(abs(wgt))),
            type="bar", orientation="h",
            y=~cluster, x=~weight, color=I("grey")) %>%
      layout(yaxis=list(showticklabels=FALSE)),
    nrows=2, heights=c(.2,.8), widths=c(.8,.2)
  ) %>% layout(showlegend=FALSE)
  return(subp)
}

