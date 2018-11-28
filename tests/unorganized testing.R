

## unorganized testing
library(swCRTdesign)
library(microbenchmark)

## DesMat
construct_DesMat(Cl=c(2,0,1))
construct_DesMat(Cl=c(2,2))

construct_DesMat(Cl=c(5,5),design="parallel")
construct_DesMat(Cl=c(5,5),design="parallel",timepoints=1)


## CovBlk
timepoints=3; sigma=2; tau=rep(2,3)
construct_CovBlk(timepoints=5, sigma=2, tau=2)
construct_CovBlk(timepoints=3, sigma=2, tau=rep(2,3))


## CovMat
construct_CovMat(SumCl=2,timepoints=3, sigma=list(c(1,2,2),c(1,1,2)),tau=0.3)
construct_CovMat(SumCl=2,timepoints=3, sigma=list(c(1,2,2),c(1,1,2)),
                 tau=list(c(.2,.2,.1),c(.2,.1,.1)))
construct_CovMat(SumCl=2,timepoints=3, sigma=list(c(1,2,2),c(1,1,2)),
                 tau=list(c(.2,.1,.1),c(.2,.2,.1)),N=c(4,4))

construct_CovMat(SumCl=2,timepoints=3, sigma=list(c(1,2,2),c(1,1,2)),
                 tau=0,N=c(1,1))
construct_CovMat(SumCl=2,timepoints=3, sigma=list(c(1,2,2),c(1,1,2)),
                 tau=0,N=c(25,16))


## wlsMixedPower
wlsMixedPower(EffSize = .1,sigma=.1,tau=.01,Cl=c(1,1,1,1))

Cl <- c(2,2,2,2,2); sigma <- 4; tau <- 0.9 ; EffSize <- 1.5
wlsMixedPower(Cl=Cl,sigma=sigma,tau=tau,EffSize=EffSize)
swPwr(swDsn(Cl),distn="gaussian",1,0,EffSize,
      tau=tau,eta=0,rho=0,sigma=sigma)


wlsMixedPower(1,Cl=c(2,0,2,0,1),1,0.5)
wlsMixedPower(EffSize=1,Cl=c(1,1,1,1,1),sigma=2 ,        tau=0.2, N=c(1,1,1,1,1) )
wlsMixedPower(EffSize=1,Cl=c(1,1,1,1,1),sigma=2*sqrt(2) ,tau=0.2, N=c(2,2,2,2,2) )

##wlsGlmmPower
wlsGlmmPower(Cl=c(1,1,1),mu0=0.04,mu1=0.02,tau=0.0,N=250)
swPwr(swDsn(c(1,1,1)),mu0=.04,mu1=.02,tau=.0,eta=0,n=250,distn="binomial")

## binomial <-> gaussian analogy
## noch benoetigt fuer den "Wrapper"
swPwr(swDsn(rep(25,4)),"binomial",200,0.03,0.025,tau=0.01,eta=0)
sigtmp <- sqrt(.0275*.9725/200)
microbenchmark(
  swPwr(swDsn(rep(25,4)),"gaussian",1,0.03,0.025,tau=0.01,eta=0,sigma=sigtmp)
  ,
  wlsMixedPower(EffSize=0.005,Cl=c(25,25,25,25),sigma=sigtmp,tau=0.01)[[1]]
  ,times=100)


## KidSafe Setting Cl
swPwr(swDsn(c(2,2,2,2,2,2),extra.time=1,tx.effect=.5), distn="binomial",
      n=250, mu0=0.03, mu1=0.02, tau=0.00262, eta=0.00131, rho=-.25, retDATA=FALSE)




