#'
#'  'construct_CovBlk'
#'
#' @param timepoints not implemented
#' @param sigma numeric, residual error of cluster means if no N given.
#' Else residual error on individual level
#' @param tau numeric, standard deviation of random intercepts
#'
#'
#'


timepoints <- 4; sigma <- 3; tau <- c(2,2,4,4)
construct_CovBlk <- function(timepoints,sigma,tau){
  return(diag(sigma^2,timepoints) + tau^2)
}


#'
#'  'construct_CovMat'
#'
#' @param I integer (vector), number of clusters per wave (in SWD)
#' @param timepoints not implemented
#' @param sigma numeric, residual error of cluster means if no N given.
#' Else residual error on individual level
#' @param tau numeric, standard deviation of random intercepts
#' @param family distribution , not implemented
#' @param N integer (vector), number of individuals per cluster.
#' Defaults to 'rep(1,sum(I))' if not passed.
#'
#'


construct_CovMat <- function(I,timepoints=NULL,sigma,tau,family=gaussian(),N=NULL){

  SumCluster <- sum(I)

  if(is.null(N)) {
    NVec <- rep(1,SumCluster)
  }else if(length(N)==SumCluster) {
    NVec <- N
  }else {stop('length of cluster sizes does not fit to total number of clusters')}

  Sigmas <- sigma / sqrt(NVec)

  if(is.null(timepoints))  timepoints <- length(I)+1   ## not implemented

  CovBlks <- mapply(construct_CovBlk,sigma=Sigmas,
                    MoreArgs=list(timepoints=timepoints,
                                  tau=tau),
                    SIMPLIFY = FALSE)
  return(Matrix::bdiag(CovBlks))
}

