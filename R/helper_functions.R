
#' @title Compute Power of a Wald Test
#'
#' computes the power of a scaled Wald test given a standard error,
#' an effect size, the degrees of freedom of the t-distribution and
#' a significance level.
#' Computes the exact power, see second example
#'
#' @param d numeric, raw effect
#' @param se numeric, standard error
#' @param df numeric, degrees of freedom of the t-distribution
#' @param sig.level numeric, significance level, defaults to 0.05
#'
#' @return a scalar
#' @export
#'
#' @examples tTestPwr(4,1,10) ; tTestPwr(4,1,30) ; tTestPwr(4,1,Inf)

tTestPwr <- function(d,se,df,sig.level=0.05){
  dsz <- abs(d/se)
  q   <- qt(sig.level/2,df=df)
  Pwr <- pt(dsz + q, df=df) + pt(-dsz + q, df=df)
  return(Pwr)
}

# tTestPwr uses a central t-distribution.
# As applied sometimes in the literature,
# e.g. Li, Turner, Preisser 2018 | Li, Redden 2015
# This is NOT the same as a t-Test, which uses a non-central t-distribution !!
# Hence the name of the above function should be changed to 'scaledWaldPwr'
# tTestPwr2 tries to mimic a t-test when only df are given.

tTestPwr2 <- function(d,se,df,sig.level=.05){
  d   <- abs(d/se)
  q   <- qt(sig.level/2, df=df)
  ncp <- -sqrt((df+2)/4)*d ## TODO: ?? how to calculate n using df ??
  pwr <- pt(q,  df, ncp=ncp) + pt(-q, df,ncp=ncp,lower.tail=FALSE)
  return(pwr)
}



################################################################################
## auxiliary functions for binomial outcome ####
tau_to_tauLin <- function(tau,mu){tau/logit.deriv(mu)}
logit.deriv   <- function(x) 1/(x-x^2)

muMarg_to_muCond <- function(muMarg,tauLin) {
  muMargLin <- binomial()$linkfun(muMarg)
  i <- function(x,muMargLin,tauLin){
    stats::dnorm(x,0,tauLin)/(1+exp(-x-muMargLin))}

  ifelse(tauLin<1e-5,muMarg,
         stats::integrate(i,-Inf,Inf,
                   muMargLin=muMargLin,tauLin=tauLin,
                   rel.tol=100*.Machine$double.eps)$value)
}
muCond_to_muMarg <- function(muCond,tauLin){
  uniroot(function(x) {muMarg_to_muCond(x,tauLin=tauLin)-muCond},
          c(1e-6,1-1e-6),
          tol=100*.Machine$double.eps)$root
}


################################################################################
## Alternative input options for covariance structure ####

### Transform icc and cac to random effects

icc_to_RandEff <- function(icc, cac=1, sigResid){
  if(any(c(icc,cac)<0,c(icc,cac)>1))
    stop("ICC and CAC must be between 0 and 1.")

  if (cac == 1) {
    gamma <- 0
    tau   <- sqrt(sigResid^2 * icc/(1-icc))
  }
  else {
    gamma <- sqrt(icc * sigResid^2 * (1-cac)/(1-icc))
    tau   <- sqrt(gamma^2 * cac / (1-cac))
  }
  return(list(gamma =gamma,
              tau   =tau))
}

RandEff_to_icc <- function(sigResid, tau, gamma=0){
  sigMarg <- sqrt(sigResid^2+tau^2+gamma^2)
  cac     <- tau^2 / (tau^2 + gamma^2)
  icc     <- (tau^2+gamma^2) / sigMarg^2
  return(list(icc=icc,
              cac=cac,
              sigMarg=sigMarg))
}


RandEff_to_alpha012 <- function(sigResid, tau, gamma, psi){
  SigMargSq <- (tau^2 + gamma^2 +psi^2 + sigResid^2)

  alpha0 <- (tau^2 + gamma^2) / SigMargSq
  alpha1 <- (tau^2) / SigMargSq
  alpha2 <- (tau^2 + psi^2)   / SigMargSq

  return(list(alpha0 = alpha0,
              alpha1 = alpha1,
              alpha2 = alpha2,
              SigMarg= sqrt(SigMargSq)))
}


alpha012_to_RandEff <- function(alpha012, sigResid=NULL, sigMarg=NULL){

  if(is.null(sigResid)==is.null(sigMarg))
    stop("Either `sigResid` or `sigMarg` must be declared (But not both).")

  a0 <- alpha012[[1]]
  a1 <- alpha012[[2]]
  a2 <- alpha012[[3]]

  sigMargSq <- if(is.null(sigResid))  sigMarg^2 else ( sigResid^2/(1-a0-a2+a1) )

  tau      <- sqrt( sigMargSq * a1 )
  gamma    <- sqrt( sigMargSq * (a0-a1) )
  psi      <- sqrt( sigMargSq * (a2-a1) )
  sigResid <- sqrt( sigMargSq * (1-a0-a2+a1) )

  return(list(tau      = tau,
              gamma    = gamma,
              psi      = psi,
              sigresid = sigResid))
}
