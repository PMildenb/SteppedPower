#'  'wlsGlmmPower'
#'
#'
#' @param I integer (vector), number of clusters per wave (in SWD)
#' @param mu0 average under control
#' @param mu1 average under intervention
#' @param tau numeric, standard deviation of random intercepts (soon) to be understood on the linear predictor level
#' @param sig.level numeric, significance level, defaults to 0.05
#' @param delay numeric (possibly vector), value between 0 and 1 specifing the
#' intervention effect in the first (second ... ) intervention phase *not implemented*
#'


## Convenience wrapper for non-normal outcomes, still in progress ...
## sigma is derived by mu_i*(1-mu_i) for control and intervention separately
## how to handle tau,eta,rho (delta-method?) ?

wlsGlmmPower <- function(I,mu0,mu1,tau,eta=NULL,rho=NULL,
                                family="binomial",N=NULL,sig.level=0.05,delay=NULL){

  DesMat <- construct_DesMat(I=I,delay=delay)
  trtmat <- matrix(DesMat[,1],nrow = sum(I),byrow=T)

  if(family =="binomial"){
    sigtmp  <- c(sigma0=sqrt(mu0*(1-mu0)),sigma1=sqrt(mu1*(1-mu1)) )
    EffSize <- mu1-mu0  ; OR <- (mu1*(1-mu0))/(mu0*(1-mu1))
    print(paste("The assumed odds ratio is",round(OR,4)))

    taus <- c(tau0=tau/(mu0*(1-mu0)),tau1=tau/(mu1*(1-mu1)))

    sigtmp <- mapply(,
                     t(trtmat), MoreArgs=list(sigma=sigtmp),SIMPLIFY = T)
    Sigmas <- split(sigtmp,rep(1:SumCl,each=timepoints))


    wlsMixedPower(DesMat=DesMat,EffSize=EffSize,
                  sigma=Sigmas,tau=taus,
                  family=family,N=N,sig.level=sig.level)
  }
}

split_sd <- function(trtvec,sd)(sd[1] + trtvec*(sd[2]-sd[1]))

