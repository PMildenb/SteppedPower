#'
#' 'zTestPower'
#'
#' @param d numeric, raw effect
#' @param se numeric, standard error
#' @param sig.level numeric, significance level, defaults to 0.05
#'


zTestPwr <- function(d,se,sig.level=0.05){
  dsz <- abs(d/se)
  Pwr <- pnorm(dsz+qnorm(sig.level/2)) + pnorm(-dsz+qnorm(sig.level/2))
  return(Pwr)
}



#'
#' 'wlsMixedPower'
#'
#' @param EffSize numeric, raw effect
#' @param I integer (vector), number of clusters per wave (in SWD)
#' @param sigma numeric, residual error of cluster means if no N given.
#' Else residual error on individual level
#' @param tau numeric, standard deviation of random intercepts
#' @param sig.level numeric, significance level, defaults to 0.05
#' @param delay numeric (possibly vector), value between 0 and 1 specifing the
#' intervention effect in the first (second ... ) intervention phase
#'

wlsMixedPower <- function(DesMat=NULL,EffSize,I,sigma,tau,family=gaussian(),
                          N=NULL,sig.level=0.05,delay=NULL){

  if(is.null(DesMat)){
    DesMat      <- construct_DesMat(I=I,delay=delay)
    timepoints  <- length(I)+1
  }
  trtmat      <- matrix(DesMat[,1],nrow = sum(I),byrow=T)

  CovMat      <- construct_CovMat(I=I,timepoints=timepoints,
                                  sigma=sigma,tau=tau,family=family,
                                  N=N,trtmat=trtmat)

  tmpmat <- t(DesMat) %*% solve(CovMat)
  VarMat <- solve(tmpmat %*% DesMat)
  HatMat <- VarMat %*% tmpmat
  Pwr    <- zTestPwr(d=EffSize,se=sqrt(VarMat[1,1]),sig.level=sig.level)
  return(list(Power=Pwr,HatMatrix=HatMat)
}

