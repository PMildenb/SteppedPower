#'  'construct_CovBlk'
#'
#' @param timepoints numeric (scalar), number of timeperiods, i.e. dimension of blocks in covariance matrix
#' @param sigma numeric (scalar or vector of length timepoints),
#' residual error (usually of cluster means)
#' @param tau numeric (scalar or vector of length timepoints), standard deviation of random intercepts,
#' A vector of length *timepoints* is interpreted as a variing sd over time (also used for binomial outcomes).
#'


construct_CovBlk <- function(timepoints,sigma,tau){
  if(length(sigma)>1 && length(sigma)!=timepoints)
    stop("length of vector sigma does not fit to number of timepoints")

  if (length(tau)==1)            taus <- rep(tau,timepoints)  else
  if (length(tau)==timepoints)   taus <- tau                  else
    stop("length of vector tau does not fit to number of timepoints")

  return(diag(sigma^2,timepoints) + taus %o% taus)
}


#'  'construct_CovMat'
#' @param SumCl total number of clusters
#' @param timepoints numeric, scalar
#' @param sigma numeric, residual error of cluster means if no N given.
#' Else residual error on individual level
#' @param tau numeric, standard deviation of random intercepts
#' @param family distribution , not implemented
#' @param N integer (vector), number of individuals per cluster.
#' Defaults to 'rep(1,sum(I))' if not passed.
#'


construct_CovMat <- function(SumCl,timepoints=NULL,sigma,tau,family=gaussian(),N=NULL){

  siglength   <- length(sigma)

  if(is.null(N)) {               NVec <- rep(1,SumCl)
  }else if(length(N)==1){        NVec <- rep(N,SumCl)
  }else if(length(N)==SumCl) {   NVec <- N
  }else stop('length of cluster sizes does not fit to total number of clusters')

  ## sigma on cluster level depending on cluster size
  if(siglength==1){
    Sigmas <- sigma / sqrt(NVec)
  }else if(siglength==timepoints){
    sigtmp <- rep(sigma,SumCl)/sqrt(rep(NVec,each=timepoints))
    Sigmas <- split(sigtmp,rep(1:SumCl,each=timepoints))
    print("sigma is assumed to change over time, but not between clusters")
  } else if(siglength==SumCl){
      if(length(sigma[[1]])==timepoints){
        Sigmas <- Map('/',sigma,sqrt(NVec))      ## used for binomial
      }else if(length(sigma[[1]])==1){
        sigtmp <- rep(sigma/sqrt(NVec),each=timepoints)
        Sigmas <- split(sigtmp,rep(1:SumCl,each=timepoints))
      }
  } else
    stop("Cannot handle length of vector sigma")


  CovBlks <- mapply(construct_CovBlk,sigma=Sigmas,tau=tau,
                    MoreArgs=list(timepoints=timepoints),
                    SIMPLIFY = FALSE)
  return(Matrix::bdiag(CovBlks))
}

