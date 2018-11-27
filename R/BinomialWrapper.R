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
#' @return
#'
#' @examples
#' wlsGlmmPower(I=c(2,2,2),mu0=0.25,mu1=0.5,tau=0.05,N=20)


## Convenience wrapper for non-normal outcomes, still in progress ...
## sigma is derived by mu_i*(1-mu_i) for control and intervention separately


wlsGlmmPower <- function(I,mu0,mu1,tau,eta=NULL,rho=NULL,
                                family="binomial",N=NULL,sig.level=0.05,delay=NULL){

  DesMat <- construct_DesMat(I=I,delay=delay)
  trtmat <- matrix(DesMat[,1],nrow = sum(I),byrow=T)
  timepoints  <- length(I)+1

  if(family =="binomial"){
    EffSize <- mu1-mu0  ; OR <- (mu1*(1-mu0))/(mu0*(1-mu1))
    print(paste("The assumed odds ratio is",round(OR,4)))

    sigma01  <- c(sigma0=sqrt(mu0*(1-mu0)),sigma1=sqrt(mu1*(1-mu1)) )
    sigtmp   <- mapply(split_sd, t(trtmat), MoreArgs=list(sd=sigma01),SIMPLIFY=T)
    Sigmas   <- split(sigtmp,rep(1:sum(I),each=timepoints))

    tau01    <- c(tau0=tau/(mu0*(1-mu0)),tau1=tau/(mu1*(1-mu1)))
    tautmp <- mapply(split_sd, t(trtmat), MoreArgs=list(sd=tau01),SIMPLIFY=T)
    Taus   <- split(tautmp,rep(1:sum(I),each=timepoints))

    wlsMixedPower(DesMat=DesMat,EffSize=EffSize,
                  timepoints=timepoints,sigma=Sigmas,tau=Taus,
                  family=family,N=N,sig.level=sig.level)
  }
}


#'  'split_sd'
#'
#' @param trtvec a vector
#' @param sd  vector of length 2 (or 3).
#'
#' @return a formula for CoxEval()
#' @export
#'
#' @examples
#' BuildFormula(VarName="Variable",AdjName="additiveAdjustment",StratName="Stratification")


split_sd <- function(trtvec,sd){
  (sd[1] + trtvec*(sd[2]-sd[1]))
}
