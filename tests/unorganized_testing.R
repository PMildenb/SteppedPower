

## unorganized testing

#### DesMat #####

## construct_time_adjust #####
construct_timeadjust(Cl=c(1,1,1), timepoints=4, "factor")
construct_timeadjust(Cl=c(1,1,1), timepoints=4, "none"  )
construct_timeadjust(Cl=c(1,1,1), timepoints=4, "linear")
tmp <- construct_timeadjust(Cl=c(1,1,1), timepoints=6, "periodic")
matplot(tmp)
tmp <- construct_timeadjust(Cl=c(1,1,1), timepoints=12, "periodic", period=6)
matplot(tmp)

## costruct_trtMat #####
construct_trtMat(Cl=c(1,1,1), trt_delay=NULL, design="SWD",      timepoints=4)
construct_trtMat(Cl=c(1,1),   trt_delay=NULL, design="parallel", timepoints=4)
construct_trtMat(Cl=c(1,1,1), trt_delay=c(.2,.4),  design="SWD",      timepoints=4)
construct_trtMat(Cl=c(1,1),   trt_delay=c(.2,.4),  design="parallel", timepoints=4)
construct_trtMat(Cl=c(2,2),   trt_delay=NULL,  design="parallel_baseline", timepoints=3)

## costruct_DesMat #####
construct_DesMat(Cl=c(2,0,1))
Des_PB <- construct_DesMat(Cl=c(2,2))

construct_DesMat(Cl=c(5,5),design="parallel",timepoints=NULL)
construct_DesMat(Cl=c(5,5),design="parallel",timepoints=1)
Des_PB <-construct_DesMat(Cl=c(2,2),design="parallel",timepoints=3)

construct_DesMat(Cl=c(2,2),design="parallel_baseline",timepoints=NULL)
construct_DesMat(Cl=c(2,2),design="parallel_baseline",timepoints=2)
Des_PB <- construct_DesMat(Cl=c(2,2),design="parallel_baseline",timepoints=c(1,2))

construct_DesMat(Cl=c(1,1,1),trt_delay=c(.3,.7))
construct_DesMat(Cl=c(1,1),  trt_delay=c(.3,.7),design="parallel")
construct_DesMat(Cl=c(1,1),  trt_delay=c(.3,.7),design="parallel_baseline")

construct_DesMat(Cl=c(1,1,1),trt_delay=c(.3,.7),time_adjust="none")
construct_DesMat(Cl=c(1,1),trt_delay=c(.3,.7),design="parallel_baseline",time_adjust="linear")

construct_DesMat(Cl=c(1,1,1,1),time_adjust="factor")
construct_DesMat(Cl=c(1,1,1,1),time_adjust="periodic")
construct_DesMat(Cl=rep(1,6),time_adjust="periodic",period=4)

## CovBlk #####
timepoints=3; sigma=2; tau=rep(2,3)
construct_CovBlk(timepoints=5, sigma=2, tau=2)
construct_CovBlk(timepoints=3, sigma=2, tau=rep(2,3))
construct_CovBlk(timepoints=3, sigma=c(1,2,3), tau=2)

eta <- c(0,.2,.2)
construct_CovBlk(timepoints,sigma,tau,eta)

## CovMat #####
construct_CovMat(SumCl=2,timepoints=3, sigma=matrix(c(1,1,2,1,2,2),nrow=2),tau=0.3)
construct_CovMat(SumCl=2,timepoints=3, sigma=1,tau=0.3,eta=1,trtMat = upper.tri(matrix(1,nrow=2,ncol=3)))

construct_CovMat(sigma=1,tau=0.3,eta=1,trtMat = upper.tri(matrix(1,nrow=2,ncol=3)))
construct_CovMat(SumCl=2,timepoints=3, sigma=1, tau=matrix(c(.3,.3,1.3,.3,1.3,1.3),nrow=2) )

construct_CovMat(SumCl=2,timepoints=3, sigma=1, tau=.3, N=4)
construct_CovMat(SumCl=2,timepoints=3, sigma=1, tau=.3, N=c(4,4))
construct_CovMat(SumCl=2,timepoints=3, sigma=1, tau=.3, N=matrix(4,2,3))

construct_CovMat(SumCl=2,timepoints=3, sigma=1, tau=.3, N=c(4,16))
construct_CovMat(SumCl=2,timepoints=3, sigma=c(1/2,1/4), tau=.3, N=1)
construct_CovMat(SumCl=2,timepoints=3, sigma=matrix(c(.5,.25),2,3), tau=.3, N=1)


Cov_PB <- construct_CovMat(4,3,sigma=matrix(1,4,3),0.1)
Cov_PB

construct_CovMat(SumCl=10,timepoints=2,sigma=1,tau=0)
construct_CovMat(SumCl=10,timepoints=2,sigma=1,tau=1)

trtMat <- construct_trtMat(c(1,1,1),NULL,"SWD")
construct_CovMat(SumCl=3,timepoints=4,sigma=sqrt(10),tau=0)
construct_CovMat(SumCl=3,timepoints=4,sigma=sqrt(10),tau=1)
construct_CovMat(SumCl=3,timepoints=4,sigma=sqrt(10),tau=1,eta=sqrt(.1),trtMat=trtMat)

## compute_wlsPower #####

## some tests for denomDF adjustment: #####
compute_wlsPower(DesMat=construct_DesMat(Cl=c(4,2,2,4)),
    EffSize=.5,sigma=1,tau=.1,N=1,df_adjust="none",sig.level=.05,verbose=T)
compute_wlsPower(DesMat=construct_DesMat(Cl=c(4,2,2,4)),
    EffSize=.5,sigma=1,tau=.1,N=1,df_adjust="between-within",sig.level=.05,verbose=T)
compute_wlsPower(DesMat=construct_DesMat(Cl=c(4,2,2,4)),
    EffSize=.5,sigma=1,tau=.1,N=1,df_adjust="containment",sig.level=.05,verbose=T)
compute_wlsPower(DesMat=construct_DesMat(Cl=c(4,4,4,4)),
    EffSize=.5,sigma=1,tau=.1,N=1,df_adjust="between-within",sig.level=.05,verbose=T)



compute_wlsPower(DesMat=construct_DesMat(Cl=c(337,337),design="parallel",timepoints=1),
    EffSize=.25,sigma=1,tau=1,N=1,df_adjust="none",sig.level=.05,verbose=F)

Cl <- rep(10,10)
DesMat <- construct_DesMat(Cl=Cl)
DesMat_prl <- construct_DesMat(Cl=c(40,40),design="parallel",timepoints=3)


SteppedPower:::compute_wlsPower(DesMat=DesMat, EffSize=.05, sigma=1,tau=.3,N=1,
                 Power=NULL,df_adjust="none",sig.level=.05,verbose=F)
wlsMixedPower(DesMat=DesMat,EffSize=.05,sigma=1,tau=.3,verbose=F)



SteppedPower:::compute_wlsPower(DesMat=DesMat_prl, EffSize=.5, sigma=1,tau=.3,N=1,
                 Power=NULL,df_adjust="none",sig.level=.05,verbose=F)
wlsMixedPower(DesMat=DesMat_prl,EffSize=.5,sigma=1,tau=.3,verbose=F)
wlsMixedPower(Cl=c(4,4),timepoints=3,design="parallel",
              EffSize=.5,sigma=1,tau=.3,verbose=F)



## optFunction #####

DesMat <- construct_DesMat(rep(1,8))
wlsMixedPower(DesMat=DesMat,EffSize=.05,sigma=1,tau=.3,N=720,verbose=F)
wlsMixedPower(DesMat=DesMat,EffSize=.05,sigma=1,tau=.3,Power=.9,verbose=F)

SteppedPower:::optFunction(DesMat=DesMat,EffSize=0.5,sigma=1,tau=.3,N=5,
                           eta=NULL,rho=NULL,Power=.9,df_adjust="none",sig.level=.05)

uniroot(SteppedPower:::optFunction,DesMat=DesMat_prl,EffSize=.5,
        sigma=1,tau=.15,Power=.9,df_adjust="none",sig.level=.05,
        lower=0.5,upper=1000)
wlsMixedPower(DesMat=DesMat_prl,EffSize=.15,sigma=1,tau=.15,N=14,verbose=F)
wlsMixedPower(DesMat=DesMat_prl,EffSize=.15,sigma=1,tau=.15,Power=.9,verbose=F,N_range = c(1,20))


## wlsMixedPower #####
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

microbenchmark::microbenchmark(
wls_parBl1 <- wlsMixedPower(EffSize=0.1,sigma=1,tau=1,Cl=c(50,50),
                           design="parallel_baseline",timepoints=51)
,
wls_parBl2 <- wlsMixedPower(EffSize=0.1,sigma=1,tau=1,Cl=c(50,50),
                           design="parallel_baseline",timepoints=c(15,36))
,
wls_swd    <- wlsMixedPower(EffSize=0.1,sigma=1,tau=1,Cl=rep(2,50),
                            design="SWD")
,times=10)
system.time(sw <- swCRTdesign::swPwr(swCRTdesign::swDsn(rep(2,50)),distn="gaussian",
                               n=1,mu0=0,mu1=.1,
                               sigma=1,tau=1,eta=0,rho=0,gamma=0))
wls_parBl1$Power
wls_parBl2$Power
wls_swd$Power
sw

wlsMixedPower(EffSize = .1,sigma=1,tau=.3,Cl=c(2,2,2,2,2),Power=.9,verbose=F)
wlsMixedPower(EffSize = .1,sigma=1,tau=.3,Cl=c(2,2,2,2,2),N=224,verbose=F)
swCRTdesign::swPwr(swCRTdesign::swDsn(c(2,2,2,2,2)),distn="gaussian",n=224,
                   mu0=0,mu1=.1,tau=.3,eta=0,rho=0,gamma=0,sigma=1)

wlsMixedPower(EffSize = .05,sigma=1,tau=.3,Cl=rep(2,20),Power=.9,verbose=F)
wlsMixedPower(EffSize = .05,sigma=1,tau=.3,Cl=rep(2,20),N=60,verbose=F)

wlsMixedPower(EffSize = .02,sigma=1,tau=.0,
              Cl=c(50,50),timepoints=5,design="parallel",
              Power=.9,verbose=F)
wlsMixedPower(EffSize = .02,sigma=1,tau=.0,
              Cl=c(10,10),timepoints=5,design="parallel",
              N=1000,verbose=F)


wlsMixedPower(Cl=c(1,1,1),trt_delay=c(.3,.7),time_adjust="None",EffSize=.1,sigma=1,tau=.1,verbose=T)
wlsMixedPower(Cl=c(1,1,1,1),trt_delay=c(.3,.7),time_adjust="None",  EffSize=.1,sigma=1,tau=.1,Power=.8,verbose=T)
wlsMixedPower(Cl=c(1,1,1,1),trt_delay=c(.3,.7),time_adjust="factor",EffSize=.1,sigma=1,tau=.01,Power=.8,verbose=T)

wlsMixedPower(Cl=c(2,2,2,2,2,2),trt_delay=c(.3,.7),time_adjust="None",EffSize=.1,sigma=1,tau=.01,Power=.8,verbose=T)
wlsMixedPower(Cl=c(2,2,2,2,2,2),trt_delay=c(.3,.7),time_adjust="factor",EffSize=.1,sigma=1,tau=.01,Power=.8,verbose=T)


wlsMixedPower(Cl=c(2,2,2,0,2,2,2,0),EffSize=.01,
              sigma=sqrt(.025*.975), time_adjust="none",
              tau=0.00254,trt_delay=.5, N=58, verbose=TRUE)
wlsMixedPower(Cl=c(2,2,2,0,2,2,2,0),EffSize=.01,
              sigma=sqrt(.025*.975), time_adjust="periodic", period=4,
              tau=0.00254,trt_delay=.5, N=58, verbose=TRUE)


## plot.wlsPower #####
plot(wlsMixedPower(Cl=c(2,5),sigma=1,tau=0.1,EffSize=1,
                            timepoints=5,design="parallel",verbose=T))


## compare_designs #####
compare_designs(EffSize=1, sigma=1 ,tau=.3, Cl=c(2,2,2,2))
compare_designs(EffSize=1, sigma=1 ,tau=.5, Cl=c(2,2,2,2))
compare_designs(EffSize=1, sigma=1 ,tau=.7, Cl=c(2,2,2,2))
compare_designs(EffSize=1, sigma=1 ,tau=1 , Cl=c(2,2,2,2))




## wlsGlmmPower #####
wlsGlmmPower(Cl=c(1,1,1),mu0=0.04,mu1=0.02,tau=0.0,N=250)
swPwr(swDsn(c(1,1,1)),mu0=.04,mu1=.02,tau=.0,eta=0,n=250,distn="binomial")

wlsGlmmPower(Cl=rep(10,5),mu0=0.04,mu1=0.02,tau=0.01,N=1, verbose=F)
swPwr(swDsn(rep(10,5)),mu0=.04,mu1=.02,tau=.01,eta=0,n=1,distn="binomial")


## binomial <-> gaussian analogy
## noch benoetigt fuer den "Wrapper"
swPwr(swDsn(rep(25,4)),"binomial",200,0.03,0.025,tau=0.01,eta=0)
sigtmp <- sqrt(.0275*.9725/200)
microbenchmark::microbenchmark(
  swPwr(swDsn(rep(25,4)),"gaussian",1,0.03,0.025,tau=0.01,eta=0,sigma=sigtmp)
  ,
  wlsMixedPower(EffSize=0.005,Cl=c(25,25,25,25),sigma=sigtmp,tau=0.01)[[1]]
  ,times=100)



## really large designs

Cl_swd <- rep(5,60) ; EffSize <- .1 ; sigma <- 1 ; tau <- 1 ; design <- "SWD"
Cl_prl <- c(150,150) ; timepoints <- 61
system.time(
  pwrSWD <- wlsMixedPower(Cl=Cl_swd, design="SWD", EffSize=EffSize, sigma=sigma, tau=tau, verbose=TRUE))
system.time(
  pwrPRL <- wlsMixedPower(Cl=Cl_prl, design="parallel", EffSize=EffSize, sigma=sigma, tau=tau,
                timepoints=timepoints, verbose=T))
dim(pwrPRL$DesignMatrix)
dim(pwrSWD$DesignMatrix)


system.time(
  HuHuPwrSWD <- swCRTdesign::swPwr(design=swCRTdesign::swDsn(rep(10,30)),distn="gaussian",
                                   n=1,mu0=0,mu1=.1,sigma=1,tau=1,eta=0,rho=0,gamma=0))

system.time(wlsMixedPower(Cl=Cl_swd, design="SWD", EffSize=EffSize, sigma=sigma, tau=tau, verbose=TRUE))


