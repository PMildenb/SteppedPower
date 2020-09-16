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
#' If design is swd, timepoints defaults to length(Cl)+1.
#' Defaults to 1 for parallel designs.
#' @param DesMat matrix of dimension , if supplied, `timepoints`,`Cl`,`trtDelay` are ignored.
#' @param trtDelay numeric (possibly vector), value(s) between 0 and 1 specifying
#' the intervention effect in the first (second ... ) intervention phase
#' @param incomplete integer, either a vector (only for SWD) or a matrix.
#' A vector defines the number of periods before and after the switch from
#' control to intervention that are observed. A matrix consists of 1's for
#' observed clusterperiods and 0's for unobserved clusterperiods.
#' @param timeAdjust character, specifies adjustment for time periods.
#' One of the following: "factor", "linear", "none", "periodic".
#' Defaults to "factor".
#' @param dsntype character, defines the type of design. Options are "SWD",
#' "parallel" and "parallel_baseline", defaults to "SWD".
#' @param mu0 numeric (scalar), mean under control
#' @param mu1 numeric (scalar), mean under treatment
#' @param marginal_mu logical. Only relevant for non-gaussian outcome.
#' Indicates whether mu0 and mu1 are to be interpreted as marginal prevalence
#' under control  and under treatment, respectively, or whether they denote
#' the prevalence conditional on random effects being 0 (It defaults to the latter).
#' @param sigma numeric, residual error of cluster means if no N given.
#' @param tau numeric, standard deviation of random intercepts
#' @param eta numeric, standard deviation of random slopes
#' @param tauAR numeric (scalar), value between 0 and 1. Defaults to NULL.
#' If `tauAR` is not NULL, the random intercept `tau` is AR1-correlated. *Currently not compatible with `rho`!=0 !*
#' @param rho numeric (scalar), correlation of `tau` and `eta`
#' @param gamma numeric (scalar), random time effect
#' @param N numeric, number of individuals per cluster. Either a scalar, vector
#' of length #Clusters or a matrix of dimension #Clusters x timepoints
#' @param family character, distribution family. One of "gaussian", "binomial".
#' Defaults to "gaussian"
#' @param Power numeric, a specified target power. If supplied, the minimal `N` is returned.
#' @param N_range numeric, vector specifying the lower and upper bound for `N`, ignored if `Power` is NULL.
#' @param sig.level numeric (scalar), significance level, defaults to 0.05
#' @param dfAdjust character, one of the following: "none","between-within", "containment", "residual".
#' @param verbose logical, should the function return the design and covariance matrix?
#' @param period numeric (scalar)
#' @param CovMat numeric, a positive-semidefinite matrix of dimension
#' (#Clusters \eqn{\cdot} timepoints) *#Cluster* \eqn{\cdot} *timepoints*
#' rows/columns. If `CovMat` is given, `sigma`, `tau`, `eta` and `rho` are ignored.
#'
#' @details see vignette 'Getting Started'
#'
#' @return an element of class wlsPower, a list that contains power and the
#' parameters of the specific setting.
#' If requested (by verbose=TRUE) it contains also relevant matrices.
#' @export
#'
#' @examples
#' ## See also vignette for more examples
#' ##
#' ##
#' ## stepped wedge design with 5 Clusters in 5 waves, residual sd = 2,
#' ## cluster effect sd = 0.2, and 10 Individuals per cluster
#' wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=10)
#' ##
#' ##
#' ## ... with auto-regressive cluster effect
#' wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, tauAR=0.7, N=10)
#' ##
#' ##
#' ## ... with varying cluster size
#' wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=c(12,8,10,9,14))
#' wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33,
#'               N=matrix(c(12,8,10,9,14,
#'                          11,8,10,9,13,
#'                          11,7,11,8,12,
#'                          10,7,10,8,11,
#'                           9,7, 9,7,11,
#'                           9,6, 8,7,11),5,6))
#' ##
#' ##
#' ## ... with random treatment effect (with sd=0.2), which is correlated with
#' ## the cluster effect with rho=0.25
#' wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, eta=.2, rho=.25, N=10)
#' ##
#' ##
#' ## ... with missing observations (a.k.a. incomplete stepped wedge design)
#' wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=10, incomplete=3)
#' wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=10,
#'              incomplete=matrix(c(1,1,1,0,0,
#'                                  1,1,1,1,0,
#'                                  1,1,1,1,1,
#'                                  1,1,1,1,1,
#'                                  0,1,1,1,1,
#'                                  0,0,1,1,1),5,6))
#'##
#'##
#'##
#'## longitudinal parallel design, with 5 time periods, 3 clusters in treatment
#'## and control arm each.
#' wlsMixedPower(mu0=0, mu1=1, Cl=c(3,3), sigma=2, tau=0.33, N=10,
#'               dsntype="parallel", timepoints=5)
#'##
#'##
#'## ... with one baseline period and four parallel periods
#' wlsMixedPower(mu0=0, mu1=1, Cl=c(3,3), sigma=2, tau=0.33, N=10,
#'               dsntype="parallel_baseline", timepoints=c(1,4))
#'##
#'##
#'##
#'## cross-over design with two timepoints before and two after the switch
#' wlsMixedPower(mu0=0, mu1=1, Cl=c(3,3), sigma=2, tau=0.33, N=10,
#'               dsntype="crossover", timepoints=c(2,2))

wlsMixedPower <- function(Cl            =NULL,
                          timepoints    =NULL,
                          DesMat        =NULL,
                          trtDelay      =NULL,
                          incomplete    =NULL,
                          timeAdjust    ="factor",
                          period        =NULL,
                          dsntype        ="SWD",
                          mu0,
                          mu1,
                          marginal_mu   =FALSE,
                          sigma         =1, ## default needed for CovMat input, still experimental
                          tau           =0,
                          eta           =NULL,
                          tauAR         =NULL,
                          rho           =NULL,
                          gamma         =NULL,
                          CovMat        =NULL,
                          N             =NULL,
                          Power         =NULL,
                          family        ="gaussian",
                          N_range       =c(1,1000),
                          sig.level     =0.05,
                          dfAdjust     ="none",
                          verbose       =FALSE){
  ## CHECKS #####
  if(!is.null(N) & !is.null(Power))
    stop("Both target power and individuals per cluster not NULL.")
  if(!is.null(rho)){
    if(is.null(eta))
      stop("If the correlation rho between random intercept and slope is 0,",
           "a random slope must be provided.")
    if( (-1)>rho | rho>1 )
      stop("Correlation rho must be between -1 and 1")
  }

  ## check DesMat #####
  if(is.null(DesMat)){
    if(!is.null(timepoints))
      if(length(trtDelay)>timepoints)
        stop("The length of vector trtDelay must be less or equal to timepoints.")

    DesMat    <- construct_DesMat(Cl         =Cl,
                                  trtDelay   =trtDelay,
                                  dsntype    =dsntype,
                                  timepoints =timepoints,
                                  timeAdjust =timeAdjust,
                                  period     =period)
  }else{
    if(min(sapply(list(Cl, timepoints, trtDelay, period),is.null))==0)
      warning("If argument DesMat is provided, Cl, timepoints, trtDelay,",
              "timeAdjust, period and dsntype are ignored.")

    if(inherits(DesMat,"matrix") & !inherits(DesMat,"DesMat")){
      DesMat <- construct_DesMat(trtmatrix=DesMat)
    }else if(!inherits(DesMat,"DesMat"))
      stop("In wlsMixedPower: Cannot interpret input for DesMat. ",
           "It must be either an object of class DesMat or a matrix")
  }

  ## incomplete designs #####
  if(!is.null(incomplete) & is.null(CovMat)){
    timepoints <- DesMat$timepoints
    lenCl      <- length(DesMat$Cl)
    SumCl      <- sum(DesMat$Cl)

    if(is.vector(incomplete) & dsntype=="SWD"){
      if(incomplete>timepoints) {
        incomplete <- timepoints
        warning("Argument `incomplete` must be less or equal to the number of",
                "timepoints. `incomplete` is set to ", timepoints )
      }
      Toep <- toeplitz(c(rep(1,incomplete),rep(Inf,lenCl-incomplete)))
      lastCols <- (timepoints-lenCl+1):timepoints

      IM <- matrix(1,lenCl,timepoints)
      IM[lower.tri(IM)]                       <- Toep[lower.tri(Toep)]
      IM[,lastCols][upper.tri(IM[,lastCols])] <- Toep[upper.tri(Toep)]

      IM <- IM[rep(1:lenCl,DesMat$Cl),]

    }else if(is.matrix(incomplete)){
      if(!nrow(incomplete) %in% c(lenCl,SumCl) | ncol(incomplete)!=timepoints)
        stop("matrix dimensions of argument `incomplete` are ",dim(incomplete),
             " but must be ", dim(DesMat$trtMat), "or ",
             dim(unique(DesMat$trtMat)))
      IM <- incomplete
      IM[which(IM==0)] <- Inf
      if(nrow(incomplete)==lenCl) IM <- IM[rep(1:lenCl,DesMat$Cl),]
    }
    sigma <- matrix(sigma, nrow=SumCl, ncol=timepoints,
                    byrow=ifelse(length(sigma)!=timepoints,TRUE,FALSE)) * IM
  }

  if(family =="binomial"){
    SumCl      <- sum(DesMat$Cl)
    timepoints <- dim(DesMat$trtMat)[2]

    MISSING_CONTROL_VARIABLE <- FALSE
    if(MISSING_CONTROL_VARIABLE){
      ## TODO: Decide whether needed or not
      ## TODO: Adjust eta, rho
      ## TODO: tau output is matrix, next if-clause "marginal_mu" demands scalar

      tau0 <- tau_to_tauLin(tau,mu0)
      tau1 <- tau_to_tauLin(tau,mu1)
      tau  <- matrix(tau0, nrow=SumCl, ncol=timepoints) + DesMat$trtMat * (tau1-tau0)
    }

    if(marginal_mu){

      mu0 <-muCond_to_muMarg(muCond=mu0, tauLin=tau)
      mu1 <-muCond_to_muMarg(muCond=mu1, tauLin=tau)
      print(paste("mu0=",round(mu0,5),", mu1=",round(mu1,5),"."))

    }

    sig0  <- sqrt(mu0*(1-mu0))
    sig1  <- sqrt(mu1*(1-mu1))
    ## for delayed trt effect only approximate sigma
    sigma <- matrix(sig0, nrow=SumCl, ncol=timepoints) + DesMat$trtMat * (sig1-sig0)

    OR <- (mu1*(1-mu0))/(mu0*(1-mu1))
    print(paste("The assumed odds ratio is",round(OR,4))) ## user information
  }

  EffSize <- mu1-mu0
  if(marginal_mu) print(paste("The (raw) effect is",EffSize))

  ## calculate samplesize (if needed) #####
  if(!is.null(Power)){
    if(Power<0 | Power>1) stop("Power needs to be between 0 and 1.")
    N_opt <- tryCatch(ceiling(
              uniroot(function(N){Power - compute_wlsPower(DesMat    =DesMat,
                                                           EffSize   =EffSize,
                                                           sigma     =sigma,
                                                           tau       =tau,
                                                           eta       =eta,
                                                           tauAR     =tauAR,
                                                           rho       =rho,
                                                           gamma     =gamma,
                                                           N         =N,
                                                           dfAdjust  =dfAdjust,
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
                          tauAR     =tauAR,
                          rho       =rho,
                          gamma     =gamma,
                          N         =N,
                          dfAdjust =dfAdjust,
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
#' @param eta numeric, standard deviation of random slopes
#' @param tauAR numeric (scalar), value between 0 and 1. Defaults to NULL. If `tauAR` is not NULL, the random intercept
#' `tau` is AR1-correlated. *Currently not compatible with `rho`!=0 !*
#' @param etaAR numeric (scalar), value between 0 and 1. Defaults to NULL. If `etaAR` is not NULL, the random slope
#' `eta` is AR1-correlated. *Currently not compatible with `rho`!=0 !*
#' @param rho numeric, correlation of tau and eta **not implemented**
#' @param gamma numeric (scalar), random time effect
#' @param N integer, number of individuals per cluster.
#' @param dfAdjust character, one of the following: **not implemented**
#' @param sig.level numeric, significance level, defaults to 0.05
#' @param verbose logical, should the function return the design and covariance matrix?
#' @param CovMat numeric, a positive-semidefinite matrix  with
#' *#Cluster* \eqn{\cdot} *timepoints* rows/columns. If `CovMat` is given, `sigma`,
#' `tau`, `eta` and `rho` are ignored.
#'
#' @export

compute_wlsPower <- function(DesMat,
                             EffSize,
                             sigma,
                             tau        =0,
                             eta        =NULL,
                             tauAR      =NULL,
                             etaAR      =NULL,
                             rho        =NULL,
                             N          =NULL,
                             CovMat     =NULL,
                             dfAdjust   ="none",
                             sig.level  =.05,
                             verbose    =FALSE){
  dsnmatrix  <- DesMat$dsnmatrix
  timepoints <- DesMat$timepoints
  SumCl      <- sum(DesMat$Cl)
  trtMat     <- DesMat$trtMat

  ## Checks ####
  if(!is.null(CovMat) & sum(sapply(c(sigma, tau, eta, rho, gamma, N),is.null))>0)
    warning("If argument CovMat is provided, sigma, tau, eta, rho, gamma and N",
            "are ignored.")

  ## get covariance matrix #####
  if(is.null(CovMat))
    CovMat   <- construct_CovMat(SumCl      =SumCl,
                                 timepoints =timepoints,
                                 sigma      =sigma,
                                 tau        =tau,
                                 eta        =eta,
                                 tauAR      =tauAR,
                                 etaAR      =etaAR,
                                 rho        =rho,
                                 trtMat     =trtMat,
                                 N          =N)

  ## matrices for power calculation #####
  tmpmat <- t(dsnmatrix) %*% Matrix::chol2inv(Matrix::chol(CovMat))
  VarMat <- Matrix::solve(tmpmat %*% dsnmatrix)
  if(verbose) ProjMat <- matrix((VarMat %*% tmpmat)[1,], nrow=SumCl, byrow=TRUE)

  ## ddf for power calculation #####
  df <- switch(dfAdjust,
               "none"           = Inf,
               "between-within" = SumCl - rankMatrix(dsnmatrix),
               "containment"    = dim(dsnmatrix)[1] - SumCl,
               "residual"       = dim(dsnmatrix)[1] - rankMatrix(dsnmatrix))
  if(df<3){
    warning(dfAdjust,"-method not applicable. No DDF adjustment used.")
    df <- Inf }

  Pwr <- tTestPwr(d=EffSize, se=sqrt(VarMat[1,1]), df=df, sig.level=sig.level)
  out <- list(Power  =Pwr,
              Params =list(N         =N,
                           sigma     =sigma, ## NOT compatible with CovMat-Input (!)
                           tau       =tau,
                           eta       =eta,
                           tauAR     =tauAR,
                           etaAR     =etaAR,
                           rho       =rho,
                           denomDF   =df,
                           dfAdjust  =dfAdjust,
                           sig.level =sig.level))
  if(verbose)
    out <- append(out,
                  list(ProjMatrix       =ProjMat,
                       DesignMatrix     =DesMat,
                       CovarianceMatrix =CovMat))
  class(out) <- append(class(out),"wlsPower")
  return(out)
}

#' print.wlsPower
#'
#' @param x object of class wlsPower
#' @param ... Arguments to be passed to methods
#'
#' @method print wlsPower
#'
#' @export
#'
#'
print.wlsPower <- function(x, ...){
  cat("Power                                = ", x$Power,    "\n")
  if(x$Params$dfAdjust!="none"){
  cat("ddf adjustment                       = ", x$Params$dfAdjust,"\n")
  cat("Denominator degrees of freedom       = ", x$Params$denomDF,  "\n")
  }
  cat("Significance level (two sided)       = ", x$Params$sig.level,"\n")

  if("N_opt" %in% names(x))
  cat("Needed N per cluster per period      = ", x$N_opt,"\n" )
}



#' plot.wlsPower
#'
#' @param x object of class wlsPower
#' @param ... Arguments to be passed to methods
#'
#' @method plot wlsPower
#'
#' @export
#'
plot.wlsPower <- function(x,...){
  if(!"ProjMatrix" %in% names(x))
    stop("Please rerun wlsMixedPower() with `verbose=TRUE` ")
  wgt <- x$ProjMatrix
  mx <- max(abs(wgt))
  sumCl <- dim(wgt)[1]
  timep <- dim(wgt)[2]
  subp <- suppressWarnings(subplot(
    plot_ly(data=data.frame(time=1:dim(wgt)[2], weight=colSums(abs(wgt))),  ## arithmetic or harmonic mean/sum ??
            type="bar", x=~time, y=~weight, color=I("grey")) %>%
      layout(yaxis=list(title=TeX('\\Sigma\\text{|weights|}')),
             xaxis=list(title="", showticklabels=FALSE))
    ,
    plotly_empty(type="scatter",mode="marker")
    ,
    plot_ly(x=1:timep,y=1:sumCl,z=wgt,type="heatmap",
            colors=grDevices::colorRamp(c("steelblue","white","firebrick")),
            xgap=.3,ygap=.3) %>%
      colorbar(len=1,limits=c(-mx,mx)) %>%
      layout(xaxis=list(title="time"),
             yaxis=list(title="cluster", autorange="reversed"))
    ,
    plot_ly(data=data.frame(cluster=1:dim(wgt)[1], weight=rowSums(abs(wgt))),
            type="bar", orientation="h",
            y=~cluster, x=~weight, color=I("grey")) %>%
      layout(xaxis=list(title=TeX('\\Sigma\\text{|weights|}')),
             yaxis=list(title="", showticklabels=FALSE, autorange="reversed"))
    ,
    nrows=2, heights=c(.2,.8), widths=c(.8,.2), titleX=TRUE, titleY=TRUE
  ) %>% layout(showlegend=FALSE) %>% config(mathjax = 'cdn'))
  return(subp)
}

