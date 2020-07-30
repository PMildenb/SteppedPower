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
#' @param N numeric, number of individuals per cluster. Either a scalar, vector
#' of length #Clusters or a matrix of dimension #Clusters x timepoints
#' @param family character, distribution family. Defaults to "gaussian" **not implemented**
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
                          time_adjust="factor", period=NULL,
                          design="SWD",
                          EffSize,
                          sigma,
                          tau,
                          eta=NULL,
                          rho=NULL,
                          CovBlk=NULL,
                          N=NULL,
                          Power=NULL,
                          family="gaussian",
                          N_range=c(1,1000),
                          sig.level=0.05,
                          df_adjust="none",
                          verbose=FALSE){

  if(!is.null(N) & !is.null(Power))
    stop("Both target power and individuals per cluster not NULL.")

  if(is.null(DesMat)){
    DesMat    <- construct_DesMat(Cl=Cl,trt_delay=trt_delay,design=design,
                                  timepoints=timepoints,time_adjust=time_adjust,
                                  period=period)
  }else
  if(inherits(DesMat,"matrix") & !inherits(DesMat,"DesMat")){
    DesMat <- construct_DesMat(trtmatrix=DesMat) ## TODO : dimension checks
  }else
  if(!inherits(DesMat,"DesMat"))
    stop("In wlsMixedPower: Cannot interpret input for DesMat. ",
         "It must be either an object of class DsnMat or a matrix")


  if(is.null(Power)){
    out <- compute_wlsPower(DesMat=DesMat,
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
                                      sigma=1,  ## It works, but why did i set this to 1 ??
                                      tau=tau,
                                      Power=Power,
                                      df_adjust=df_adjust,
                                      sig.level=.05,
                                      interval=N_range)$root),
                      error=function(cond){
                                     message(paste0("Maximal N yields power below ",Power,
                                                    ". Increase argument N_range."))
                                     return(N_range[2])})
    out <- compute_wlsPower(DesMat=DesMat,
                            EffSize=EffSize,
                            sigma=sigma,
                            tau=tau,
                            N=N_opt,
                            sig.level=sig.level,
                            df_adjust=df_adjust,
                            CovBlk=CovBlk,
                            verbose=verbose)
    out$N_opt <- N_opt
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

  diff <- (Power - compute_wlsPower(DesMat=DesMat,
                                    EffSize=EffSize,
                                    sigma=sigma,
                                    tau=tau,
                                    N=N,
                                    df_adjust=df_adjust,
                                    sig.level=sig.level,
                                    CovBlk=NULL,
                                    verbose=FALSE)$Power)
  return(diff)}



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
#'
#' @export

compute_wlsPower <- function(DesMat,
                             EffSize,
                             sigma,
                             tau,
                             N,
                             Power,
                             CovBlk=NULL,
                             df_adjust=df_adjust,
                             sig.level,
                             verbose){
  dsnmatrix  <- DesMat$dsnmatrix
  timepoints <- DesMat$timepoints
  SumCl      <- DesMat$SumCl


  if(is.null(N)) N <- 1
  NMat <- if(length(N) %in% c(1,timepoints,SumCl*timepoints)) {
            matrix(N, nrow=SumCl, ncol=timepoints)
          }else stop(paste('length of cluster size vector N is ', N,
                     '. This does not fit to given number of clusters, which is ',SumCl))

  lenS <- length(sigma)
  sigmaMat <- if(lenS==1){                      matrix(sigma, nrow=SumCl, ncol=timepoints)
              }else if(lenS==timepoints){       matrix(sigma, nrow=SumCl, ncol=timepoints)
              }else if(lenS==SumCl){            matrix(sigma, nrow=SumCl, ncol=timepoints, byrow=TRUE)
              }else if(lenS==timepoints*SumCl) {matrix(sigma, nrow=SumCl, ncol=timepoints)
              }else stop(paste('length of sigma is ', N,
                               '. This does not fit to given number of timepoints, which is ',timepoints,
                               ' or to the given number of Clusters, which is ', SumCl))
  if(timepoints==SumCl & lenS==SumCl) warning("sigma is assumed to change over time. If you wanted sigma
                                              to change between clusters, please provide as matrix of dimension
                                              #Cluster x timepoints")

  sigmaMat <- sigmaMat / sqrt(NMat)
  sigmaLst <- split(sigmaMat,row(sigmaMat))

  # if(!is.null(rho)){
  #   matrix(dsnmatrix$matrix[,"trtvec"],
  # }

  CovMat   <- construct_CovMat(timepoints=timepoints, sigma=sigmaLst, tau=tau)

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

  out <- list(Power     =tTestPwr(d=EffSize, se=sqrt(VarMat[1,1]), df=df, sig.level=sig.level),
              denomDF   =df,
              df_adjust =df_adjust,
              sig.level =sig.level)
  if(verbose) out <- append(out,list(
                WeightMatrix=WgtMat,
                DesignMatrix=dsnmatrix,
                CovarianceMatrix=CovMat))
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
  cat("ddf adjustment                       = ", x$df_adjust,"\n")
  if(x$df_adjust!="none")
  cat("Denominator degrees of freedom       = ", x$denomDF,  "\n")

  cat("Significance level (two sided)       = ", x$sig.level,"\n")

  if("N_opt" %in% names(x))
  cat("Needed N per cluster per period      = ", x$N_opt,"\n" )
}



#' plot.wlsPower
#'
#' @param x w
#' @param ... Arguments to be passed to methods
#'
#' @method plot wlsPower
#'
#' @export
#'
plot.wlsPower <- function(x, ...){

  WgtMat <- x$WeightMatrix
  HatData <- reshape2::melt(t(WgtMat[c(nrow(WgtMat):1),]))
  names(HatData) <- c("Time","Cluster","Weight")

  plotraw <- ggplot2::ggplot(HatData,ggplot2::aes_string("Time","Cluster")) +
    ggplot2::theme_minimal()  + ggplot2::scale_fill_gradient2(low="steelblue",mid="white",high="red") +
    ggplot2::geom_tile(ggplot2::aes_string(fill="Weight"),colour="white")


  Timeweights    <- colSums(abs(WgtMat))
  plotCluster <- ggplot2::ggplot(data.frame(Cluster=1:dim(WgtMat)[1],
                                            Clusterweights=rowSums(abs(WgtMat))),
                                 ggplot2::aes_string("Cluster","Clusterweights")) +

    ggplot2::geom_point() + ggplot2::theme_minimal()
  plotPeriods <- ggplot2::ggplot(data.frame(Periods=1:dim(WgtMat)[2],
                                            Weights=colSums(abs(WgtMat))),
                                 ggplot2::aes_string("Periods","Weights")) +
    ggplot2::geom_point() + ggplot2::theme_minimal()

  return(list(plotraw,plotCluster,plotPeriods))
}

