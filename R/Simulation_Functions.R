## Simulation functions


#' SimSlim
#'
#' @param design as in swDsn()
#' @param n Number per Cluster
#' @param mu0 mean response under control
#' @param mu1 mean response under treatment
#' @param time.effect a vector of fixed time (i.e. period) effects
#' @param sigma residual variance on indiviual level
#' @param tau standard deviation of random cluster effect
#' @param eta standard deviation of random treatment effect
#' @param rho correlation of rho tau and eta
#' @param time.lab defaults to NULL
#' @param seed a fixed seed
#' @param seed.designonly logical, if TRUE, only random cluster/treatment effect is determined by the seed.
#'  If false, realisations of response as well.
#'
#' @export


## simuliert fuer gegebenes SWD binomiale outcomes (cluster*time aggregiert)
Simslim <- function (design, n, mu0, mu1, time.effect, sigma, tau,
                     eta, rho = 0, time.lab = NULL, seed = NULL, seed.designonly= F)
{
  distn <- "binomial"
  link <- "logit"
  p0 <- mu0
  p1 <- mu1
  theta <- p1 - p0
  CV.p0 <- tau/p0
  CV.theta <- eta/abs(theta)
  X.ij <- design$swDsn
  if (length(time.effect) == 1) {
    time.effectVec <- rep(time.effect, design$total.time)
  } else
    if (length(time.effect) == design$total.time) {
      time.effectVec <- time.effect
    } else {
      warning("Invalid time.effects length. Either specify a scalar fixed time
              effect (i.e., the same fixed time effect at each time point),
              or specify a vector of fixed time effects for each time point.
              It is best to ignore any results which follow as a result of this
              warning message, and correctly assign value(s) to the
              'time.effect' function argument of swSim().")
    }
  beta.ij <- matrix(rep(time.effectVec, design$n.clusters),
                    design$n.clusters, design$total.time, byrow = TRUE)
  thetaX.ij <- X.ij * theta
  sigMat <- matrix(c(tau^2, rho * tau * eta, rho * tau * eta,
                     eta^2), 2, 2)
  set.seed(seed)
  if (tau != 0 & eta != 0) {
    z <- rnorm(2 * design$n.clusters)
    zMat <- matrix(z, nrow = 2, byrow = FALSE)
    ar.i <- chol(sigMat) %*% zMat
    a.i <- rep(ar.i[1, ], design$total.time)
    r.i <- rep(ar.i[2, ], design$total.time)
  } else
    if (tau != 0 & eta == 0) {
      a.i <- rep(rnorm(design$n.clusters, 0, tau), design$total.time)
      r.i <- 0
    } else
      if (tau == 0 & eta != 0) {
        a.i <- 0
        r.i <- rep(rnorm(design$n.clusters, 0, eta), design$total.time)
      } else
        if (tau == 0 & eta == 0) {
          a.i <- 0
          r.i <- 0
        }
  a.ij <- matrix(a.i, nrow = design$n.clusters, ncol = design$total.time,
                 byrow = FALSE)
  rX.ij <- X.ij * matrix(r.i, nrow = design$n.clusters, ncol = design$total.time,
                         byrow = FALSE)


  link.mu.ij <- mu0 + beta.ij + thetaX.ij + a.ij + rX.ij

  if (is.null(time.lab)) {
    time.lab <- 1:design$total.time
    time.ij <- matrix(rep(c(time.lab), design$n.clusters),
                      design$n.clusters, design$total.time, byrow = TRUE)
  } else
    if (length(time.lab) == design$total.time) {
      time.ij <- matrix(rep(c(time.lab), design$n.clusters),
                        design$n.clusters, design$total.time, byrow = TRUE)
    } else
    {
      warning("The length of the specified 'time.lab' vector is not equal to the total time points determined based on the SW design from specifying 'cluters', 'extra.time', and 'all.ctl.time0'. Ignore any results that follow and re-specify the vector 'time.lab' if you want to specify the label at each time point; specify 'NULL' to have the default labeling 1, 2, ... to the total number of time points.")
    }

  cluster.ij <- matrix(rep(c(1:design$n.clusters), each = design$total.time),
                       design$n.clusters, design$total.time, byrow = TRUE)

  mu.vec <- 1/(1 + exp(-(as.vector(t(link.mu.ij))))) ## inv.logit
  if(seed.designonly) set.seed(NULL)
  Y.ij <- rbinom(design$total.time * design$n.clusters , n , mu.vec )

  response.var    <- Y.ij
  tx.var          <- as.vector(t(X.ij))
  time.var        <- as.factor(as.vector(t(time.ij)))       ## 'as.factor' added
  cluster.var     <- as.factor(as.vector(t(cluster.ij)))

  swData      <- data.frame(response.var, tx.var, time.var, cluster.var)

  return(swData)
}



## GENERAL SIMULATION FORMULA (cluster effect on linear predictor)

#' realpower.simulate
#'
#' @param n helper
#'
#'
#' @export

realpower.simulate <- function(n=1,nInd=100,tau,eta=0,rho=0,mu0,mu1,design,
                               whichModel=c("AGQ","PQL","GEE"),nAGQ=1){

  sizes      <- rep(nInd,design$n.clusters*design$total.time)
  SimHHshort <- Simslim(design,n=nInd,car::logit(mu0),car::logit(mu1),
                        tau=tau,eta=eta,rho=rho,time.effect=0)
  out.p      <- NULL
  out.tx     <- NULL

  if("AGQ" %in% whichModel){
    modAGQ <- lme4::glmer(response.var/nInd ~ tx.var  + time.var + (1|cluster.var),
                   weights =sizes,
                   family = "binomial" ,data=SimHHshort)
    p.ztest.AGQ <- coef(summary(modAGQ))["tx.var","Pr(>|z|)"]
    tx.est.AGQ  <- modAGQ@beta[2]

    out.p  <- c(out.p,p.ztest.AGQ=p.ztest.AGQ)
    out.tx <- c(out.tx,tx.est.AGQ=tx.est.AGQ)
  }

  if("PQL" %in% whichModel){
    modPQL <- MASS::glmmPQL(response.var/sizes ~ tx.var + time.var,
                           random=~1|cluster.var,
                           weights = sizes,
                           family = "binomial",data=SimHHshort,verbose = FALSE)
    p.ttest.PQL <- coef(summary(modPQL))["tx.var","p-value"]
    tx.est.PQL  <- as.numeric(modPQL$coefficients$fixed[2])
    denDF       <- coef(summary(modPQL))["tx.var","DF"]

    out.p  <- c(out.p,p.ttest.PQL=p.ttest.PQL)
    out.tx <- c(out.tx,tx.est.PQL=tx.est.PQL)
  }

  if("GEE" %in% whichModel){
    modGEE <- suppressMessages(gee(
      cbind(response.var,nInd-response.var) ~ tx.var + time.var,
      id = cluster.var,
      # corstr = "exchangeable",
      corstr = "AR-M", Mv=1,
      maxiter = 1000,
      family = "binomial",SimHHshort))
    tx.est.GEE   <- as.numeric(modGEE$coefficients[2])
    test         <- coef(summary(modGEE))["tx.var",c("Estimate","Robust S.E.")]
    p.ztest.GEE  <- as.numeric(2*pnorm(-abs(test[1]),0,test[2]))

    out.p  <- c(out.p,p.ztest.GEE=p.ztest.GEE)
    out.tx <- c(out.tx,tx.est.GEE=tx.est.GEE)
  }

  ifelse("PQL" %in% whichModel,
         out <- c(out.p,out.tx,denDF=denDF),
         out <- c(out.p,out.tx))
  return( out )
}
