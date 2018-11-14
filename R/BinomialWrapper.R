#'  'binomialPower'
#'
#'
#'
#'
#'
#'
#'
#'
#'

## Convenience wrapper for non-normal outcomes, in progress ...


wlsGlmmPower <- function(I,mu0,mu1,tau,eta=NULL,rho=NULL,
                                family=binomial(),N=NULL,sig.level=0.05){

  Sigmas  <- c(sigma0=mu0*(1-mu0),sigma1=mu1*(1-mu1))
  EffSize <- mu1-mu0


  wlsMixedPower(EffSize=EffSize,I=I,sigma=Sigmas,tau=tau,
                family=family,N=N,sig.level=sig.level)
}


wlsGlmmPower(I=c(1,1,1),mu0=0.03,mu1=0.02,tau=0.015)
