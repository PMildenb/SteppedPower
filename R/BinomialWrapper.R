#'  wlsGlmmPower
#'
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD)
#' @param mu0 average under control
#' @param mu1 average under intervention
#' @param tau numeric, standard deviation of random intercepts (soon) to be understood on the linear predictor level
#' @param eta *not implemented*
#' @param rho *not implemented*
#' @param design study design. One of the following values: "SWD", "parallel", "parallel_baseline" ...
#' @param family currently only "binomial"
#' @param N Number of Individuals per cluster/period
#' @param Power Target power. If between 0 and 1, function returns needed `N`
#' @param sig.level numeric, significance level, defaults to 0.05
#' @param delay numeric (possibly vector), value between 0 and 1 specifing the
#' intervention effect in the first (second ... ) intervention phase *not implemented*
#' @param verbose
#'
#' @return as `wlsMixedPower`
#' @export
#'
#' @examples
#' wlsGlmmPower(Cl=c(2,2,2),mu0=0.25,mu1=0.5,tau=0.05,N=20)


## Convenience wrapper for non-normal outcomes, still in progress ...
## sigma is derived by mu_i*(1-mu_i) for control and intervention separately


wlsGlmmPower <- function(Cl,mu0,mu1,tau,eta=NULL,rho=NULL,design="SWD",timepoints=NULL,
                         family="binomial",N=NULL,sig.level=0.05,delay=NULL,
                         verbose=TRUE){

  if(is.null(timepoints)){
    if(design=="SWD"){      timepoints  <- length(Cl)+1 } else
    if(design=="parallel"){ timepoints  <- 1            } else
    if(design=="parallel_baseline") {timepoints <- 2    }
  }
  DesMat <- construct_DesMat(Cl=Cl,delay=delay,design=design,timepoints=timepoints)
  trtmat <- matrix(DesMat[,1],nrow = sum(Cl),byrow=T)

  if(family =="binomial"){
    EffSize <- mu1-mu0  ; OR <- (mu1*(1-mu0))/(mu0*(1-mu1))
    print(paste("The assumed odds ratio is",round(OR,4)))

    sigma01  <- c(sigma0=sqrt(mu0*(1-mu0)),sigma1=sqrt(mu1*(1-mu1)) )
    sigtmp   <- mapply(SteppedPower::split_sd, t(trtmat), MoreArgs=list(sd=sigma01),SIMPLIFY=T)
    Sigmas   <- split(sigtmp,rep(1:sum(Cl),each=timepoints))

#    tau01    <- c(tau0=tau/(mu0*(1-mu0)),tau1=tau/(mu1*(1-mu1)))
    tau01  <- c(tau0=tau,tau1=tau)
    tautmp <- mapply(SteppedPower::split_sd, t(trtmat), MoreArgs=list(sd=tau01),SIMPLIFY=T)
    Taus   <- split(tautmp,rep(1:sum(Cl),each=timepoints))

    wlsMixedPower(DesMat=DesMat,EffSize=EffSize,
                  timepoints=timepoints,sigma=Sigmas,tau=Taus,
                  family=family,N=N,Power=NULL,sig.level=sig.level,verbose=verbose)
  }
}


#'  split_sd
#'
#'  small helper function for wlsGlmmPower. returns sd_0 for periods in control and sd_1
#'  for interventional periods for each cluster.
#'
#' @param trtvec a vector
#' @param sd  vector of length 2 (or 3).
#'
#' @return numeric, scalar.
#' @export
#'

split_sd <- function(trtvec,sd){
  (sd[1] + trtvec*(sd[2]-sd[1]))
}
