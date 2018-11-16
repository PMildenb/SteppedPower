#'  'wlsGlmmPower'
#'
#'
#' @param I integer (vector), number of clusters per wave (in SWD)
#' @param
#' @param tau numeric, standard deviation of random intercepts (soon) to be understood on the linear predictor level
#' @param sig.level numeric, significance level, defaults to 0.05
#' @param delay numeric (possibly vector), value between 0 and 1 specifing the
#' intervention effect in the first (second ... ) intervention phase *not implemented
#'
#'
#'
#'
#'
#'

## Convenience wrapper for non-normal outcomes, in progress ...
## how to handle tau,eta,rho (delta-method?) ?

wlsGlmmPower <- function(I,mu0,mu1,tau,eta=NULL,rho=NULL,
                                family="binomial",N=NULL,sig.level=0.05,delay=NULL){
  if(family =="binomial"){
    Sigmas  <- c(sigma0=sqrt(mu0*(1-mu0)),sigma1=sqrt(mu1*(1-mu1)) )
    EffSize <- mu1-mu0  ; OR <- (mu1*(1-mu0))/(mu0*(1-mu1))
    print(paste("The assumed odds ratio is",round(OR,4)))

    wlsMixedPower(EffSize=EffSize,I=I,sigma=Sigmas,tau=tau,
                family=family,N=N,sig.level=sig.level,delay=delay)
  }
}


wlsGlmmPower(I=c(2,2,2,2,2,2,0),mu0=0.03,mu1=0.02,tau=0.0,N=250)



## binomial <-> gaussian analogy
## noch benoetigt fuer den "Wrapper"
swPwr(swDsn(rep(25,4)),"binomial",200,0.03,0.025,tau=0.01,eta=0)
sigtmp <- sqrt(.0275*.9725/200)
microbenchmark(
  swPwr(swDsn(rep(25,4)),"gaussian",1,0.03,0.025,tau=0.01,eta=0,sigma=sigtmp)
  ,
  wlsMixedPower(EffSize=0.005,I=c(25,25,25,25),sigma=sigtmp,tau=0.01)
  ,times=100)
