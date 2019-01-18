## Simulation functions

library(car)
library(gee)
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

realpower.simulate <- function(n=1,nInd=100,tau,eta=0,rho=0,mu0,mu1,design,
                               whichModel=c("HQ","PQL","GEE")){

  sizes      <- rep(nInd,design$n.clusters*design$total.time)
  SimHHshort <- Simslim(design,n=nInd,car::logit(mu0),car::logit(mu1),
                        tau=tau,eta=eta,rho=rho,time.effect=0)
  out.p      <- NULL
  out.tx     <- NULL

  if("HQ" %in% whichModel){
    mod1a <- glmer(response.var/nInd ~ tx.var  + time.var + (1|cluster.var),
                   weights =sizes,
                   family = "binomial" ,data=SimHHshort)
    p.ztest.1a <- coef(summary(mod1a))["tx.var","Pr(>|z|)"]
    tx.est.1a  <- mod1a@beta[2]

    out.p  <- c(out.p,p.ztest.1a=p.ztest.1a)
    out.tx <- c(out.tx,tx.est.1a=tx.est.1a)
  }

  if("PQL" %in% whichModel){
    mod3a <- MASS::glmmPQL(response.var/sizes ~ tx.var + time.var,
                           random=~1|cluster.var,
                           weights = sizes,
                           family = "binomial",data=SimHHshort,verbose = FALSE)
    p.ttest.3a <- coef(summary(mod3a))["tx.var","p-value"]
    tx.est.3a  <- as.numeric(mod3a$coefficients$fixed[2])
    denDF      <- coef(summary(mod3a))["tx.var","DF"]

    out.p  <- c(out.p,p.ttest.3a=p.ttest.3a)
    out.tx <- c(out.tx,tx.est.3a=tx.est.3a)
  }

  if("GEE" %in% whichModel){
    mod5a <- invisible(gee(
      cbind(response.var,nInd-response.var) ~ tx.var + time.var,
      id = cluster.var,
      corstr = "exchangeable",
      #     corstr = "unstructured",
      maxiter = 1000,
      family = "binomial",SimHHshort))
    tx.est.5a   <- as.numeric(mod5a$coefficients[2])
    test      <- coef(summary(mod5a))["tx.var",c("Estimate","Robust S.E.")]
    p.ztest.5a  <- as.numeric(2*pnorm(-abs(test[1]),0,test[2]))

    out.p  <- c(out.p,p.ztest.5a=p.ztest.5a)
    out.tx <- c(out.tx,tx.est.5a=tx.est.5a)
  }

  ifelse("PQL" %in% whichModel,
         out <- c(out.p,out.tx,denDF=denDF),
         out <- c(out.p,out.tx))
  return( out )
}
