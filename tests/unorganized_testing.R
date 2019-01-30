

## unorganized testing
library(swCRTdesign)
library(microbenchmark)




#### DesMat

## construct_time_adjust
construct_timeadjust(Cl=c(1,1,1), timepoints=4, "factor")
construct_timeadjust(Cl=c(1,1,1), timepoints=4, "none"  )
construct_timeadjust(Cl=c(1,1,1), timepoints=4, "linear")

## costruct_trtvec
construct_trtvec(Cl=c(1,1,1), trt_delay=NULL, design="SWD",      timepoints=4)
construct_trtvec(Cl=c(1,1),   trt_delay=NULL, design="parallel", timepoints=4)
construct_trtvec(Cl=c(1,1,1), trt_delay=c(.2,.4),  design="SWD",      timepoints=4)
construct_trtvec(Cl=c(1,1),   trt_delay=c(.2,.4),  design="parallel", timepoints=4)

## costruct_DesMat
construct_DesMat(Cl=c(2,0,1))
Des_PB <- construct_DesMat(Cl=c(2,2))

construct_DesMat(Cl=c(5,5),design="parallel",timepoints=NULL)
construct_DesMat(Cl=c(5,5),design="parallel",timepoints=1)
Des_PB <-construct_DesMat(Cl=c(2,2),design="parallel",timepoints=3)

construct_DesMat(Cl=c(2,2),design="parallel_baseline",timepoints=NULL)
construct_DesMat(Cl=c(2,2),design="parallel_baseline",timepoints=2)
Des_PB <- construct_DesMat(Cl=c(2,2),design="parallel_baseline",timepoints=c(1,2))

construct_DesMat(Cl=c(1,1,1),trt_delay=c(.3,.7))
construct_DesMat(Cl=c(1,1),trt_delay=c(.3,.7),design="parallel")
construct_DesMat(Cl=c(1,1),trt_delay=c(.3,.7),design="parallel_baseline")


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
Cov_PB <- construct_CovMat(4,3,1,0.1)


## wlsInnerFunction
Cl <- rep(10,10)
DesMat <- construct_DesMat(Cl=Cl)
DesMat_prl <- construct_DesMat(Cl=c(40,40),design="parallel",timepoints=3)


wlsInnerFunction(DesMat=DesMat, EffSize=.05,
                 sigma=1,tau=.3,N=1,
                 Power=NULL,sig.level=.05,verbose=F)
wlsMixedPower(DesMat=DesMat,EffSize=.05,sigma=1,tau=.3,verbose=F)


wlsInnerFunction(DesMat=DesMat_prl, EffSize=.5,
                 sigma=1,tau=.3,N=1,
                 Power=NULL,sig.level=.05,verbose=F)
wlsMixedPower(DesMat=DesMat_prl,EffSize=.5,sigma=1,tau=.3,verbose=F)
wlsMixedPower(Cl=c(4,4),timepoints=3,design="parallel",
              EffSize=.5,sigma=1,tau=.3,verbose=F)


## optFunction

wlsMixedPower(DesMat=DesMat,EffSize=.05,sigma=1,tau=.3,N=46,verbose=F)
wlsMixedPower(DesMat=DesMat,EffSize=.05,sigma=1,tau=.3,Power=.9,verbose=F)

SteppedPower:::optFunction(DesMat=DesMat_prl,EffSize=0.5,
            sigma=1,tau=.3,N=1,
            Power=.9,sig.level=.05)

uniroot(SteppedPower:::optFunction,DesMat=DesMat_prl,EffSize=.5,
        sigma=1,tau=.15,Power=.9,sig.level=.05,
        lower=0.5,upper=1000)
wlsMixedPower(DesMat=DesMat_prl,EffSize=.5,sigma=1,tau=.15,N=14,verbose=F)
wlsMixedPower(DesMat=DesMat_prl,EffSize=.5,sigma=1,tau=.15,Power=.9,verbose=F)


## wlsMixedPower
wlsMixedPower(EffSize = .1,sigma=1,tau=.01,Cl=c(1,1,1,1),verbose=F)

Cl <- c(2,2,2,2,2); sigma <- 4; tau <- 0.9 ; EffSize <- 1.5
wlsMixedPower(Cl=Cl,sigma=sigma,tau=tau,EffSize=EffSize,verbose=F)
swCRTdesign::swPwr(swCRTdesign::swDsn(Cl),distn="gaussian",1,0,EffSize,
      tau=tau,eta=0,rho=0,sigma=sigma)

wlsMixedPower(EffSize=1,Cl=c(1,1,1,1,1),sigma=2 ,        tau=0.2, N=c(1,1,1,1,1) )
wlsMixedPower(EffSize=1,Cl=c(1,1,1,1,1),sigma=2*sqrt(2) ,tau=0.2, N=c(2,2,2,2,2) )

Cl=c(3,3) ; EffSize=1 ; sigma=1 ; tau=0.1 ; family=gaussian() ; timepoints=1 ;
design="parallel" ; N <- NULL ; sig.level=0.05
wlsMixedPower(Cl=c(5,5),EffSize=.2,sigma=1,tau=0.1,timepoints=2,design="parallel")


wls_parBl1 <- wlsMixedPower(EffSize=0.1,sigma=1,tau=1,Cl=c(50,50),
                           design="parallel_baseline",timepoints=51)
wls_parBl2 <- wlsMixedPower(EffSize=0.1,sigma=1,tau=1,Cl=c(50,50),
                           design="parallel_baseline",timepoints=c(15,36))
wls_swd    <- wlsMixedPower(EffSize=0.1,sigma=1,tau=1,Cl=rep(2,50),
                            design="SWD")
wls_parBl1$Power
wls_parBl2$Power
wls_swd$Power

wlsMixedPower(EffSize = .1,sigma=1,tau=.3,Cl=c(2,2,2,2,2),Power=.9,verbose=F)
wlsMixedPower(EffSize = .1,sigma=1,tau=.3,Cl=c(2,2,2,2,2),N=224,verbose=F)
swPwr(swDsn(c(2,2,2,2,2)),distn="gaussian",n=224,mu0=0,mu1=.1,tau=.3,eta=0,sigma=1)

wlsMixedPower(EffSize = .05,sigma=1,tau=.3,Cl=rep(2,20),Power=.9,verbose=F)
wlsMixedPower(EffSize = .05,sigma=1,tau=.3,Cl=rep(2,20),N=60,verbose=F)

wlsMixedPower(EffSize = .02,sigma=1,tau=.0,
              Cl=c(50,50),timepoints=5,design="parallel",
              Power=.9,verbose=F)
wlsMixedPower(EffSize = .02,sigma=1,tau=.0,
              Cl=c(10,10),timepoints=5,design="parallel",
              N=1000,verbose=F)


wlsMixedPower(Cl=c(1,1,1),trt_delay=c(.3,.7),time_adjust="None",EffSize=.1,sigma=1,tau=.1,verbose=T)
wlsMixedPower(Cl=c(1,1,1),trt_delay=c(.3,.7),time_adjust="None",  EffSize=.1,sigma=1,tau=.1,Power=.8,verbose=T)
wlsMixedPower(Cl=c(1,1,1,1),trt_delay=c(.3,.7),time_adjust="factor",EffSize=.1,sigma=1,tau=.01,Power=.8,verbose=T)

wlsMixedPower(Cl=c(2,2,2,2,2,2),trt_delay=c(.3,.7),time_adjust="None",EffSize=.1,sigma=1,tau=.01,Power=.8,verbose=T)
wlsMixedPower(Cl=c(2,2,2,2,2,2),trt_delay=c(.3,.7),time_adjust="factor",EffSize=.1,sigma=1,tau=.01,Power=.8,verbose=T)


## wlsGlmmPower
wlsGlmmPower(Cl=c(1,1,1),mu0=0.04,mu1=0.02,tau=0.0,N=250)
swPwr(swDsn(c(1,1,1)),mu0=.04,mu1=.02,tau=.0,eta=0,n=250,distn="binomial")


wlsGlmmPower(Cl=rep(10,5),mu0=0.04,mu1=0.02,tau=0.01,N=1, verbose=F)
swPwr(swDsn(rep(10,5)),mu0=.04,mu1=.02,tau=.01,eta=0,n=1,distn="binomial")


## plot_wlsPower
plot_wlsPower(wlsMixedPower(1,Cl=c(2,5),1,0.1,family="gaussian",
                            timepoints=NULL,design="parallel",verbose=T))
plot_wlsPower(wls_parBl2)[[1]]



## compare_designs
compare_designs(EffSize=1, sigma=1 ,tau=.3, Cl=c(2,2,2,2))
compare_designs(EffSize=1, sigma=1 ,tau=.5, Cl=c(2,2,2,2))
compare_designs(EffSize=1, sigma=1 ,tau=.7, Cl=c(2,2,2,2))
compare_designs(EffSize=1, sigma=1 ,tau=1 , Cl=c(2,2,2,2))




## binomial <-> gaussian analogy
## noch benoetigt fuer den "Wrapper"
swPwr(swDsn(rep(25,4)),"binomial",200,0.03,0.025,tau=0.01,eta=0)
sigtmp <- sqrt(.0275*.9725/200)
microbenchmark(
  swPwr(swDsn(rep(25,4)),"gaussian",1,0.03,0.025,tau=0.01,eta=0,sigma=sigtmp)
  ,
  wlsMixedPower(EffSize=0.005,Cl=c(25,25,25,25),sigma=sigtmp,tau=0.01)[[1]]
  ,times=100)






## try switch
N <- rep(2,6)
SumCl <- 6

Cl <- c(2,2) ; timepoints <- 3 ; time_adjust <- "none"



## really large designs

Cl_swd <- rep(10,60) ; EffSize <- .01 ; sigma <- 1 ; tau <- 0 ; design <- "SWD"
Cl_prl <- c(300,300) ; timepoints <- 61
system.time(
  wlsMixedPower(Cl=Cl_swd, design="SWD", EffSize=EffSize, sigma=sigma, tau=tau))
system.time(
  wlsMixedPower(Cl=Cl_prl, design="SWD", EffSize=EffSize, sigma=sigma, tau=tau,
                timepoints=timepoints, verbose=T))

