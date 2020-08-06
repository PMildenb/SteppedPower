
#' construct_CovBlk
#'
#' constructs the covariance matrix for multiple measurements of one cluster.
#' This function is not designed to be used directly.
#'
#' @param timepoints numeric (scalar), number of timeperiods, Cl.e. dimension of blocks in covariance matrix
#' @param sigma numeric (scalar or vector of length `timepoints`),
#' residual error (usually of cluster means)
#' @param tau numeric (scalar or vector of length *timepoints*), standard deviation of random intercepts,
#' A vector of length *timepoints* is interpreted as a variing sd over time (also used for binomial outcomes).
#' @param eta numeric (vector of length *timepoints*), standard deviation of random slope
#' @param rho numeric (scalar), correlation of random effects `tau` and `eta`.
#'
#' @return a block of a covariance matrix, corresponding to one cluster
#' @export
#'
#' @examples
#' construct_CovBlk(timepoints=5, sigma=2, tau=2)


construct_CovBlk <- function(timepoints,
                             sigma,
                             tau,
                             eta=NULL,
                             rho=NULL){
  # if(length(sigma)>1 && length(sigma)!=timepoints)
  #   stop("length of vector sigma does not fit to number of timepoints")
  if (length(tau)==1)            taus <- rep(tau,timepoints)  else
  if (length(tau)==timepoints)   taus <- tau                  else
    stop("length of vector tau does not fit to number of timepoints")

  out <- diag(sigma^2,timepoints) + taus %o% taus
  if(!is.null(eta)) out <- out + eta %o% eta
  if(!is.null(rho)) {
    rhoVec <- (rho * eta * tau)
    out    <- out + outer(rhoVec, rhoVec, "+")
  }
  return(out)
}


#' construct_CovMat
#'
#' constructs a covariance matrix, a bloc diagonal matrix. Calls `construct_CovBlk`
#' for each block.
#'
#' @param SumCl total number of clusters
#' @param timepoints numeric, scalar
#' @param sigma numeric (scalar, vector or matrix), residual error of cluster means.
#' @param tau numeric (scalar, vector or matrix), standard deviation of random intercepts
#' @param N numeric (scalar, vector or matrix), number of individuals per cluster.
#' Defaults to 'rep(1,sum(Cl))' if not passed.
#' @param eta numeric (scalar), standard deviation of random slopes. If `eta` is
#' given, `trtMat` is needed as well.
#' @param rho numeric (scalar), correlation of random effects `tau` and `eta`.
#' @param trtMat a matrix of dimension *#Cluster* x *timepoints* as produced by
#' the function `construct_trtMat`, indicating the cluster-periods that receive
#' interventional treatment. Defaults to NULL. If trtMat is given, the arguments
#' `SumCl` and `timepoints` are ignored (!).
#' @param CovBlk a matrix of dimension *timepoints* x *timepoints*.
#'
#' @return a covariance matrix
#' @export
#'
#' @examples
#' construct_CovMat(SumCl=2,timepoints=3,
#'                  sigma=list(c(1,2,2),c(1,1,2)),
#'                  tau=list(c(.2,.1,.1),c(.2,.2,.1)),N=c(20,16))




construct_CovMat <- function(SumCl      =NULL,
                             timepoints =NULL,
                             sigma,
                             tau,
                             eta        =NULL,
                             rho        =NULL,
                             gamma      =NULL,
                             trtMat     =NULL,
                             N          =NULL,
                             CovBlk     =NULL){
  if(!is.null(CovBlk)){
    CovBlks <- rep(list(CovBlk),SumCl)
  } else {
    ## Checks ##
    if(is.null(SumCl)      & !is.null(trtMat)) SumCl      <- nrow(trtMat)
    if(is.null(timepoints) & !is.null(trtMat)) timepoints <- ncol(trtMat)
    timepoints  <- sum(timepoints)  ## dirty hack to fix potential vector input for parll+baseline

    if(!is.null(eta) & is.null(trtMat)) stop("If `eta` is not NULL, a treatment matrix must be given to construct_CovMat()")
    if(!is.null(rho) & is.null(eta))    stop("In construct_CovMat: eta is needed if rho is not NULL")

    ## N ##
    if(is.null(N)) N <- 1
    NMat <- if(length(N) %in% c(1,SumCl,SumCl*timepoints)) {
      matrix(N, nrow=SumCl, ncol=timepoints)
    }else stop(paste('length of cluster size vector N is ', N,
                     '. This does not fit to given number of clusters, which is ',SumCl))

    ## sigma ##
    lenS <- length(sigma)
    sigmaMat <- if(lenS==1){          matrix(sigma, nrow=SumCl, ncol=timepoints)
    }else if(lenS==SumCl){            matrix(sigma, nrow=SumCl, ncol=timepoints)
    }else if(lenS==timepoints){       matrix(sigma, nrow=SumCl, ncol=timepoints, byrow=TRUE)
    }else if(lenS==timepoints*SumCl) {matrix(sigma, nrow=SumCl, ncol=timepoints)
    }else stop(paste('length of sigma is ', N,
                     '. This does not fit to given number of timepoints, which is ',timepoints,
                     ' or to the given number of clusters, which is ', SumCl))
    if(timepoints==SumCl & lenS==SumCl) warning("sigma is assumed to change over time. If you wanted sigma
                                              to change between clusters, please provide as matrix of dimension
                                              #Cluster x timepoints")

    ## N into sigma ##
    sigmaMat <- sigmaMat / sqrt(NMat)
    sigmaLst <- split(sigmaMat,row(sigmaMat))

    ## tau , eta and rho  ##
    tauMat <- matrix(tau, nrow=SumCl, ncol=timepoints) ## tau can be scalar, vector or matrix
    tauLst <- split(tauMat, row(tauMat))


    if(!is.null(eta)) {
      etaMat <- trtMat * eta                           ## eta input must be a scalar
      etaLst <- split(etaMat, row(etaMat))             ## eta is passed as vector of length SumCl
    } else etaLst <- vector("list", length=SumCl)

    if(!is.null(rho)) {
      rhoLst <- as.list(rep(rho,SumCl))                ## rho input must be a scalar
    } else rhoLst <- vector("list", length=SumCl)      ## rho is passed as scalar

    CovBlks <- mapply(construct_CovBlk,
                      sigma = sigmaLst,
                      tau   = tauLst,
                      eta   = etaLst,
                      rho   = rhoLst,
                      MoreArgs = list(timepoints=timepoints),
                      SIMPLIFY = FALSE)
  }
  CovMat <- Matrix::bdiag(CovBlks)
  if(!is.null(gamma)) {
    M <- matrix(1, nrow=SumCl, ncol=SumCl)
    M <- M %x% diag(timepoints)*gamma
    CovMat <- CovMat + M
  }

  return(CovMat)
}
