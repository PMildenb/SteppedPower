

## unorganized testing
library(swCRTdesign)
library(microbenchmark)

## DesMat
construct_DesMat(I=c(2,0,1))
construct_DesMat(I=c(2,2))

## CovBlk
timepoints=3; sigma=2; tau=rep(2,3)
construct_CovBlk(timepoints=5, sigma=2, tau=2)
construct_CovBlk(timepoints=3, sigma=2, tau=rep(2,3))


## CovMat
construct_CovMat(I=c(1,1), sigma=1,      tau=0.3)
construct_CovMat(I=c(1,1), sigma=c(1,2), tau=0.3,
                 trtmat=matrix(c(0,0,1,0,1,1),nrow=2))

SumCl <- 2; timepoints <- 3; sigma <- split()
CovMat      <- construct_CovMat(SumCl=SumCl,timepoints=timepoints,
                                sigma=sigma,tau=tau,family=family,
                                N=N,trtmat=trtmat)


## wlsMixedPower
wlsMixedPower(EffSize = .1,sigma=.1,tau=.01,I=c(1,1,1,1))

I <- c(2,2,2,2,2); sigma <- 4; tau <- 0.9 ; EffSize <- 1.5
wlsMixedPower(I=I,sigma=sigma,tau=tau,EffSize=EffSize)
swPwr(swDsn(I),distn="gaussian",1,0,EffSize,
      tau=tau,eta=0,rho=0,sigma=sigma)


wlsMixedPower(1,I=c(2,0,2,0,1),1,0.5)
wlsMixedPower(EffSize=1,I=c(1,1,1,1,1),sigma=2 ,        tau=0.2, N=c(1,1,1,1,1) )
wlsMixedPower(EffSize=1,I=c(1,1,1,1,1),sigma=2*sqrt(2) ,tau=0.2, N=c(2,2,2,2,2) )

##wlsGlmmPower
wlsGlmmPower(I=c(1,1,1),mu0=0.03,mu1=0.02,tau=0.0,N=250)

## binomial <-> gaussian analogy
## noch benoetigt fuer den "Wrapper"
swPwr(swDsn(rep(25,4)),"binomial",200,0.03,0.025,tau=0.01,eta=0)
sigtmp <- sqrt(.0275*.9725/200)
microbenchmark(
  swPwr(swDsn(rep(25,4)),"gaussian",1,0.03,0.025,tau=0.01,eta=0,sigma=sigtmp)
  ,
  wlsMixedPower(EffSize=0.005,I=c(25,25,25,25),sigma=sigtmp,tau=0.01)
  ,times=100)


## KidSafe Setting I
swPwr(swDsn(c(2,2,2,2,2,2),extra.time=1,tx.effect=.5), distn="binomial",
      n=250, mu0=0.03, mu1=0.02, tau=0.00262, eta=0.00131, rho=-.25, retDATA=FALSE)





##  understanding swPwr
View(swPwr)
J <- timepoints <- 3 ; tau <- 0.5 ; sigSq <- 4 ; n <- 10 ; I <- SumCluster <- 4 ; eta <- 0.04
I.rep <- waves <- c(2,2) ; i.indx <- 1 ; swDesignUnique <- matrix(c(0,1,1,0,0,1),2,byrow=T)

Wmat.blk <- tau^2 + diag(sigSq/n, J)
Wmat.partial <- kronecker(diag(1, I), Wmat.blk)
Xij.Xil.ARRAY <- array(NA, c(J, J, length(I.rep)))
for (i.indx in 1:length(I.rep)) {
  Xi.jl <- swDesignUnique[i.indx, ]
  Xij.Xil.ARRAY[, , i.indx] <- (Xi.jl %o% Xi.jl) * eta^2 +
    outer(Xi.jl, Xi.jl, "+") * rho * eta * tau
}
Xij.Xil.LIST_blk <- sapply(1:length(I.rep), function(x) NULL)
