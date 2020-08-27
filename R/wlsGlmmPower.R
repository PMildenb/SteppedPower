#'  wlsGlmmPower
#'
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD)
#' @param mu0 numeric, (marginal) average under control
#' @param mu1 numeric, (marginal) average under intervention
#' @param marginal_mu logical,
#' @param tau numeric, standard deviation of random intercepts, interpreted as normal noise on the proportions
#' @param tau_lin numeric, standard deviation of random intercepts on the linear predictor
#' @param eta *not implemented*
#' @param rho *not implemented*
#' @param design study design. One of the following values: "SWD", "parallel", "parallel_baseline" ...
#' @param family currently only "binomial"
#' @param N Number of Individuals per cluster/period
#' @param Power Target power. If between 0 and 1, function returns needed `N`
#' @param sig.level numeric, significance level, defaults to 0.05
#' @param dfAdjust character, one of the following: **not implemented**
#' @param delay numeric (possibly vector), value between 0 and 1 specifing the
#' dose of intervention in the first (second ... ) intervention phase *not implemented*
#' @param verbose logical, shall the function be chatty about what it did?
#'
#' @return as `wlsMixedPower`
#' @export
#'
# @examples
#

## Convenience wrapper for non-normal outcomes, **still experimental**
## sigma is derived by mu_i*(1-mu_i) for control and intervention separately


wlsGlmmPower <- function(Cl            =NULL,
                         timepoints    =NULL,
                         DesMat        =NULL,
                         trtDelay      =NULL,
                         incomplete    =NULL,
                         timeAdjust    ="factor",
                         period        =NULL,
                         design        ="SWD",
                         mu0,
                         mu1,
                         marginal_mu   =TRUE,
                         tau           =0,
                         tau_lin       =NULL,
                         eta           =NULL,
                         tauAR         =NULL,
                         rho           =NULL,
                         gamma         =NULL,
                         CovMat        =NULL,
                         N             =NULL,
                         Power         =NULL,
                         family        ="binomial",
                         N_range       =c(1,1000),
                         sig.level     =0.05,
                         dfAdjust     ="none",
                         verbose       =FALSE){

  ## CHECKS #####
  if(!is.null(N) & !is.null(Power))
    stop("Both target power and individuals per cluster not NULL.")
  if(!is.null(rho)){
    if(is.null(eta))
      stop("If the correlation rho between random intercept and slope is 0, a random slope must be provided.")
    if( (-1)>rho | rho>1 )
      stop("Correlation rho must be between -1 and 1")
  }
  if(!is.null(DesMat)){
    if(min(sapply(list(Cl, timepoints, trtDelay, incomplete, period),is.null))==0)
      warning("If argument DesMat is provided, Cl, timepoints, trtDelay, incomplete, timeAdjust, period and design are ignored.")
    else{
      if(!is.null(timepoints))
        if(length(trtDelay)>timepoints)
          stop("The length of vector trtDelay must be less or equal to timepoints.")
    }
  }

  ## check DesMat #####
  if(is.null(DesMat)){
    DesMat    <- construct_DesMat(Cl         =Cl,
                                  trtDelay   =trtDelay,
                                  design     =design,
                                  timepoints =timepoints,
                                  timeAdjust =timeAdjust,
                                  period     =period)
  }else if(inherits(DesMat,"matrix") & !inherits(DesMat,"DesMat")){
    DesMat <- construct_DesMat(trtmatrix=DesMat) ## TODO : dimension checks
  }else if(!inherits(DesMat,"DesMat"))
    stop("In wlsMixedPower: Cannot interpret input for DesMat. ",
         "It must be either an object of class DesMat or a matrix")

  ## incomplete designs
  if(!is.null(incomplete) & is.null(CovMat)){
    timepoints <- DesMat$timepoints
    lenCl      <- length(DesMat$Cl)

    if(is.vector(incomplete) & design=="SWD"){
      Toep <- toeplitz(c(rep(1,incomplete),rep(Inf,lenCl-incomplete)))
      lastCols <- (timepoints-lenCl+1):timepoints

      IM <- matrix(1,lenCl,timepoints)
      IM[lower.tri(IM)]                       <- Toep[lower.tri(Toep)]
      IM[,lastCols][upper.tri(IM[,lastCols])] <- Toep[upper.tri(Toep)]

      IM <- IM[rep(1:lenCl,DesMat$Cl),]
    } else if(is.matrix(incomplete)){
      if(nrow(incomplete)!=sum(DesMat$Cl) | ncol(incomplete)!=timepoints)
        stop("matrix dimensions of argument incoplete are ",dim(incomplete), " but must be ",
             dim(DesMat$trtMat))
      IM <- incomplete
      IM[which(IM==0)] <- Inf
    }
    sigma <- matrix(sigma, nrow=lenCl, ncol=timepoints,
                    byrow=ifelse(length(sigma)!=timepoints,TRUE,FALSE)) * IM
  }


  if(family =="binomial"){

    if(MISSING_CONTROL_VARIABLE){
      tau0 <- tau_to_tauLin(tau_lin,mu0)
      tau1 <- tau_to_tauLin(tau_lin,mu1)
      tau  <- matrix(tau0, nrow=SumCl, ncol=timepoints) + DesMat$trtMat * (tau1-tau0)

      ## TODO: Adjust eta
    }

    if(marginal_mu){
      mu0 <-muCond_to_muMarg(muCond=mu0, tauLin=tau0)
      mu1 <-muCond_to_muMarg(muCond=mu1, tauLin=tau1)
      print(paste("mu0_marg=",mu0,", mu1_marg=",mu1,"."))

      sig0  <- sqrt(mu0*(1-mu0))
      sig1  <- sqrt(mu1*(1-mu1))
      sigma <- matrix(sig0, nrow=SumCl, ncol=timepoints) + DesMat$trtMat * (sig1-sig0)

      EffSize <- mu1-mu0
      print(paste("The effect size is",EffSize))

    }

    EffSize <- mu1-mu0
    OR <- (mu1_marg*(1-mu0_marg))/(mu0_marg*(1-mu1_marg))
    print(paste("The assumed odds ratio is",round(OR,4))) ## user information
  }

  ## calculate Power #####
  out <- compute_wlsPower(DesMat    =DesMat,
                          EffSize   =EffSize,
                          sigma     =sigma,
                          tau       =tau,
                          eta       =eta,
                          tauAR     =tauAR,
                          rho       =rho,
                          gamma     =gamma,
                          N         =N,
                          dfAdjust =dfAdjust,
                          sig.level =sig.level,
                          CovMat    =CovMat,
                          verbose   =verbose)
  return(out)
}

## aux functions for binomial
tau_to_tauLin <- function(tau,mu){tau/logit.deriv(mu)}
logit.deriv   <- function(x) 1/(x-x^2)
logit         <- function(x) log(x/(1-x))
inv_logit     <- function(x) exp(x)/(1+exp(x))

muMarg_to_muCond <- function(muMarg,tauLin) {
  muMargLin <- logit(muMarg)
  i <- function(x,muMargLin,tauLin){
    dnorm(x,0,tauLin)/(1+exp(-x-muMargLin))}

  ifelse(tauLin<1e-5,muMarg,
         integrate(i,-Inf,Inf,
                   muMargLin=muMargLin,tauLin=tauLin,
                   rel.tol=100*.Machine$double.eps)$value)
}
muCond_to_muMarg <- function(muCond,tauLin){
  uniroot(function(x) {muMarg_to_muCond(x,tauLin=tauLin)-muCond},
          c(1e-6,1-1e-6),
          tol=100*.Machine$double.eps)$root
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
