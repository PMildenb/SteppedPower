
#' construct_CovBlk
#'
#' constructs the covariance matrix for multiple measurements of one cluster
#'
#' @param timepoints numeric (scalar), number of timeperiods, Cl.e. dimension of blocks in covariance matrix
#' @param sigma numeric (scalar or vector of length `timepoints`),
#' residual error (usually of cluster means)
#' @param tau numeric (scalar or vector of length timepoints), standard deviation of random intercepts,
#' A vector of length *timepoints* is interpreted as a variing sd over time (also used for binomial outcomes).
#'
#' @return a block of a covariance matrix
#' @export
#'
#' @examples
#' construct_CovBlk(timepoints=5, sigma=2, tau=2)


construct_CovBlk <- function(timepoints,
                             sigma,
                             tau){
  # if(length(sigma)>1 && length(sigma)!=timepoints)
  #   stop("length of vector sigma does not fit to number of timepoints")
  if (length(tau)==1)            taus <- rep(tau,timepoints)  else
  if (length(tau)==timepoints)   taus <- tau                  else
    stop("length of vector tau does not fit to number of timepoints")

  return(diag(sigma^2,timepoints) + taus %o% taus)
}


#' construct_CovMat
#'
#' constructs a covariance matrix, a bloc diagonal matrix. Calls `construct_CovBlk`
#' for each block.
#'
#' @param SumCl total number of clusters
#' @param timepoints numeric, scalar
#' @param sigma numeric, residual error of cluster means.
#' @param tau numeric, standard deviation of random intercepts
#' @param N integer (vector), number of individuals per cluster.
#' Defaults to 'rep(1,sum(Cl))' if not passed.
#'
#' @return a covariance matrix
#' @export
#'
#' @examples
#' construct_CovMat(SumCl=2,timepoints=3,
#'                  sigma=list(c(1,2,2),c(1,1,2)),
#'                  tau=list(c(.2,.1,.1),c(.2,.2,.1)),N=c(20,16))




construct_CovMat <- function(SumCl=NULL,
                             timepoints=NULL,
                             sigma,
                             tau,
                             N=NULL,
                             CovBlk=NULL){
  if(!is.null(CovBlk)){
    CovBlks <- rep(list(CovBlk),SumCl)
    return(Matrix::bdiag(CovBlks))
  } else {
  #   timepoints  <- sum(timepoints)  ## dirty hack to fix potential vector input for parll+baseline
  #   siglength   <- length(sigma)
  #
  #   if(is.null(N)) {               NVec <- rep(1,SumCl)
  #   }else if(length(N)==1){        NVec <- rep(N,SumCl)
  #   }else if(length(N)==SumCl) {   NVec <- N
  #   # }else if(length(N)==SumCl*timepoints){ NVec <- splitN
  #   }else stop('length of cluster sizes does not fit to total number of clusters')
  #
  #   ## sigma on cluster level depending on cluster size
  #   if(siglength==1){
  #     Sigmas <- sigma / sqrt(NVec)
  #   }else if(siglength==timepoints){
  #     sigtmp <- rep(sigma,SumCl)/sqrt(rep(NVec,each=timepoints))
  #     Sigmas <- split.default(sigtmp,rep(1:SumCl,each=timepoints))
  #     print("sigma is assumed to change over time, but not between clusters")
  #   }else if(siglength==SumCl){
  #     if(length(sigma[[1]])==timepoints){
  #        Sigmas <- Map('/',sigma,sqrt(NVec))      ## used for binomial
  #     }else if(length(sigma[[1]])==1){
  #        sigtmp <- rep(sigma/sqrt(NVec),each=timepoints)
  #        Sigmas <- split.default(sigtmp,rep(1:SumCl,each=timepoints))
  #     }
  #   }else
  #     stop("Cannot handle length of vector sigma")

    CovBlks <- mapply(construct_CovBlk,sigma=sigma,tau=tau,
                      MoreArgs=list(timepoints=timepoints),
                      SIMPLIFY = FALSE)
    # CovMat  <- list(Matrix=Matrix::bdiag(CovBlks))
    # class(CovMat) <- append(class(CovMat),"CovMat")
    CovMat  <- Matrix::bdiag(CovBlks)

    return(CovMat)
  }
}

