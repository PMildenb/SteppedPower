#'  wlsGlmmPower
#'
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD)
#' @param mu0 numeric, (marginal) average under control
#' @param mu1 numeric, (marginal) average under intervention
#' @param tau numeric, standard deviation of random intercepts, interpreted as normal noise on the proportions
#' @param tau_lin numeric, standard deviation of random intercepts on the linear predictor
#' @param eta *not implemented*
#' @param rho *not implemented*
#' @param design study design. One of the following values: "SWD", "parallel", "parallel_baseline" ...
#' @param family currently only "binomial"
#' @param N Number of Individuals per cluster/period
#' @param Power Target power. If between 0 and 1, function returns needed `N`
#' @param sig.level numeric, significance level, defaults to 0.05
#' @param df_adjust character, one of the following: **not implemented**
#' @param delay numeric (possibly vector), value between 0 and 1 specifing the
#' dose of intervention in the first (second ... ) intervention phase *not implemented*
#' @param verbose logical, shall the function be chatty about what it did?
#'
#' @return as `wlsMixedPower`
#' @export
#'
#' @examples
#' wlsGlmmPower(Cl=c(2,2,2),mu0=0.25,mu1=0.5,tau=0.05,N=20)


## Convenience wrapper for non-normal outcomes, **still experimental**
## sigma is derived by mu_i*(1-mu_i) for control and intervention separately


wlsGlmmPower <- function(Cl,mu0,mu1,
                         marginal_mu=TRUE,
                         tau=NULL,tau_lin=NULL,
                         eta=NULL,rho=NULL,
                         design="SWD",timepoints=NULL,family="binomial",
                         N=NULL,sig.level=0.05,
                         trt_delay=NULL,time_adjust="factor",
                         df_adjust="none",verbose=TRUE){

  if(is.null(timepoints)){
    if(design=="SWD"){      timepoints  <- length(Cl)+1 } else
    if(design=="parallel"){ timepoints  <- 1            } else
    if(design=="parallel_baseline") {timepoints <- 2    }
  }
  DesMat <- construct_DesMat(Cl=Cl,trt_delay=trt_delay,design=design,
                             timepoints=timepoints,time_adjust=time_adjust)$matrix
  trtmat <- matrix(DesMat[,1],nrow = sum(Cl),byrow=T)

  if(family =="binomial"){
    logit.deriv <- function(x){1/(x-x^2)}
    inv_logit   <- function(x) exp(x)/(1+exp(x))

    if(marginal_mu){

      pmarg<-function(mu,tau_lin) { i<-function(x,mu,tau_lin){ dnorm(x,0,tau_lin)/(1+exp(-x-mu))}
        ifelse(tau_lin<1e-3,1/(1+exp(-mu)),integrate(i,-Inf,Inf,mu=mu,tau_lin=tau_lin)$value)}
      pcond_to_pmarglin <- function(pcond,tau_lin,low=-10){
        uniroot(function(x) {pmarg(x,tau_lin=tau_lin)-pcond},c(low,0))$root}

      mu0_marg <-inv_logit(pcond_to_pmarglin(mu0,tau=tau_lin) )
      mu1_marg <-inv_logit(pcond_to_pmarglin(mu1,tau=tau_lin) )
      print(paste("mu0_marg=",mu0,", mu1_marg=",mu1,"."))
    }

    EffSize <- mu1-mu0  ; OR <- (mu1_marg*(1-mu0_marg))/(mu0_marg*(1-mu1_marg))
    print(paste("The assumed odds ratio is",round(OR,4))) ## user information

    sigma01  <- c(sigma0=sqrt(mu0*(1-mu0)),sigma1=sqrt(mu1*(1-mu1)) )
    sigtmp   <- mapply(split_sd, t(trtmat), MoreArgs=list(sd=sigma01),SIMPLIFY=T)
    Sigmas   <- split(sigtmp,rep(1:sum(Cl),each=timepoints))

    if(!is.null(tau_lin)){
      tau.lin_to_tau.expit <- function(tau_lin,mu){tau_lin/logit.deriv(mu)}
      tau01  <- c(tau0=tau.lin_to_tau.expit(tau_lin,mu0),
                  tau1=tau.lin_to_tau.expit(tau_lin,mu1))
    } else {
      tau01  <- c(tau0=tau,tau1=tau)
    }
    print(tau01)
    tautmp <- mapply(split_sd, t(trtmat), MoreArgs=list(sd=tau01),SIMPLIFY=T)
    Taus   <- split(tautmp,rep(1:sum(Cl),each=timepoints))

    wlsMixedPower(Cl=Cl, timepoints=timepoints, DesMat=DesMat, trt_delay=trt_delay,
                  design=design, EffSize=EffSize, sigma=Sigmas,tau=Taus,
                  N=N,Power=NULL,df_adjust=df_adjust,sig.level=sig.level,verbose=verbose)
  }
}


#'  split_sd
#'
#'  small auxiliary function for wlsGlmmPower. returns sd_0 for periods in control and sd_1
#'  for interventional periods for each cluster. Not to be called directly.
#'
#' @param trtvec a vector
#' @param sd  vector of length 2 (or 3).
#'
#' @return numeric, scalar.
#'

split_sd <- function(trtvec,sd){
  (sd[1] + trtvec*(sd[2]-sd[1]))
}
