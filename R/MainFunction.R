#'
#'
#'
#'
#'
#'
#'
#'
I <- c(2,2,2,2,2); sigma <- 4; tau <- 0.9 ; EffSize <- 1.5
d <- EffSize ; se <- sqrt(VarTrt)


zTestPwr <- function(d,se,sig.level=0.05){
  dsz <- abs(d/se)
  Pwr <- pnorm(dsz+qnorm(sig.level/2)) + pnorm(-dsz+qnorm(sig.level/2))
  return(Pwr)
}



wlsMixedPower <- function(EffSize,I,sigma,tau,family=gaussian(),N=1,sig.level=0.05){

  DesMat <- construct_DesMat(I=I)
  CovMat <- construct_CovMat(I=I,sigma=sigma,tau=tau,family=family,N=N)

  VarTrt <- solve(t(DesMat) %*% solve(CovMat) %*% DesMat)[1,1]

  Pwr    <- zTestPwr(d=EffSize,se=sqrt(VarTrt),sig.level=sig.level)
  return(Pwr)
}

swPwr(swDsn(c(2,2,2,2,2)),distn = "gaussian",1,0,1.5,0.9,0,0,sigma = 4)
View(swPwr)
