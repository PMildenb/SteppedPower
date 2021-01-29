
library(SteppedPower)

## unorganized testing

################################################################################
#### DesMat #####

## construct_timeAdjust #####
construct_timeadjust(Cl=c(1,1,1), timepoints=4, "factor")
construct_timeadjust(Cl=c(1,1,1), timepoints=4, "none"  )
construct_timeadjust(Cl=c(1,1,1), timepoints=4, "linear")
tmp <- construct_timeadjust(Cl=c(1,1,1), timepoints=6, "periodic")
matplot(tmp)
tmp <- construct_timeadjust(Cl=c(1,1,1), timepoints=12, "periodic", period=6)
matplot(tmp)

construct_timeadjust(3,timeBlk=matrix(c(0,0,1,0,1,1),nrow=2))

## costruct_trtMat #####
construct_trtMat(Cl=c(1,1,1), trtDelay=NULL, dsntype="SWD",      timepoints=4)
construct_trtMat(Cl=c(1,1),   trtDelay=NULL, dsntype="parallel", timepoints=4)
construct_trtMat(Cl=c(1,1,1), trtDelay=c(.2,.4),  dsntype="SWD",      timepoints=4)
construct_trtMat(Cl=c(1,1),   trtDelay=c(.2,.4),  dsntype="parallel", timepoints=4)
construct_trtMat(Cl=c(2,2),   trtDelay=NULL,  dsntype="parallel_baseline", timepoints=3)
construct_trtMat(Cl=c(2,2),   trtDelay=NULL,  dsntype="crossover", timepoints=c(2,2))
# debugonce(construct_trtMat)

adist("b", c("SWD","parallel","parallel_baseline","crossover"),
      cost=c(insertions=1, deletions=100, substitutions=100))


## costruct_DesMat #####
construct_DesMat(Cl=c(2,0,1))
Des_PB <- construct_DesMat(Cl=c(2,2))

construct_DesMat(Cl=c(5,5),dsntype="parallel",timepoints=NULL)
construct_DesMat(Cl=c(5,5),dsntype="parallel",timepoints=1)
Des_PB <-construct_DesMat(Cl=c(2,2),dsntype="parallel",timepoints=3)

construct_DesMat(Cl=c(2,2),dsntype="parallel_baseline",timepoints=NULL)
construct_DesMat(Cl=c(2,2),dsntype="parallel_baseline",timepoints=2)
Des_PB <- construct_DesMat(Cl=c(2,2),dsntype="parallel_baseline",timepoints=c(1,2))

construct_DesMat(Cl=c(1,1,1),trtDelay=c(.3,.7))
construct_DesMat(Cl=c(1,1),  trtDelay=c(.3,.7),dsntype="parallel")
construct_DesMat(Cl=c(1,1),  trtDelay=c(.3,.7),dsntype="parallel_baseline")

construct_DesMat(Cl=c(1,1,1),trtDelay=c(.3,.7),timeAdjust="none")
construct_DesMat(Cl=c(1,1),trtDelay=c(.3,.7),dsntype="parallel_baseline",timeAdjust="linear")

construct_DesMat(Cl=c(1,1,1,1),timeAdjust="factor")
construct_DesMat(Cl=c(1,1,1,1),timeAdjust="periodic")
construct_DesMat(Cl=rep(1,6),timeAdjust="periodic",period=4)

construct_DesMat(Cl=c(1,1,2),trtDelay=.5,timeBlk=diag(4))

# debugonce(construct_DesMat)
a <- construct_DesMat(Cl=c(3,2,1),N=c(1,1,1,2,2,3))
a$dsnmatrix

b <- construct_DesMat(Cl=c(1,1,1),N=c(1,2,3),
                 trtmatrix=matrix(c(0,0,0,1,0,0,1,1,0,1,1,1),3))

c <- construct_DesMat(Cl=c(1,1,1),N=c(1,2,3), INDIV_LVL = TRUE,
                 trtmatrix=matrix(c(0,0,0,1,0,0,1,1,0,1,1,1),3))

################################################################################
## CovBlk #####
construct_CovBlk(sigma=2:4, tau=rep(1,3))
construct_CovBlk(sigma=1:3, tau=c(5,1,0) )
construct_CovBlk(sigma=c(1,2,3), tau=c(2,1,2))

eta <- c(0,.2,.2)
construct_CovBlk(sigma=rep(2,3),tau=rep(1,3),eta=eta)

rho <- .3
construct_CovBlk(sigma=rep(2,3),tau=rep(1,3),eta=eta,rho=rho)

## CovSubMat #####

tp <- 4
SteppedPower:::construct_CovSubMat(N          =3,
                                   timepoints =tp,
                                   sigma      =rep(sqrt(10000),tp),
                                   tau        =rep(sqrt(5),tp),
                                   psi        =rep(sqrt(2000),tp),
                                   eta        =c(0,0,rep(sqrt(40),(tp-2))),
                                   tauAR      =NULL,
                                   gamma      =rep(sqrt(300),tp))
SteppedPower:::construct_CovSubMat(N          =3,
                                   timepoints =tp,
                                   sigma      =rep(sqrt(10000),tp),
                                   tau        =rep(sqrt(5),tp),
                                   psi        =rep(sqrt(2000),tp),
                                   eta        =c(0,0,rep(sqrt(40),(tp-2))),
                                   tauAR      =NULL,
                                   gamma      =rep(sqrt(300),tp),
                                   INDIV_LVL  =TRUE )


## CovMat #####
construct_CovMat(SumCl=2,timepoints=3, sigma=0,tau=0.3)

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
construct_CovMat(SumCl=3,timepoints=4,sigma=sqrt(10),tau=1,eta=sqrt(.1),rho=.5,trtMat=trtMat)


construct_CovMat(SumCl=3, timepoints=4, sigma=3, tau=2, tauAR=.8)

construct_CovMat(SumCl=3, timepoints=4, sigma=3, tau=2, tauAR=0)
construct_CovMat(SumCl=3, timepoints=4, sigma=3, tau=0, gamma=2)
construct_CovMat(SumCl=3, timepoints=4, sigma=3, tau=0, gamma=2)


construct_CovMat(2,3,1,10,eta= .1, trtMat=matrix(c(0,0,1,0,1,1),2,3))
construct_CovMat(2,3,1,10,eta= .1*matrix(c(0,0,1,0,1,1),2,3))

## psi ##

# debugonce(construct_CovMat)
# debugonce(construct_CovSubMat)

construct_CovMat(SumCl=2,timepoints=4,sigma=1,tau=0,N=c(2,3),psi=5,gamma=9)
construct_CovMat(SumCl=2,timepoints=4,sigma=1,tau=0,N=c(2,3),psi=5,gamma=9,INDIV_LVL=TRUE)
construct_CovMat(SumCl=2,timepoints=4,
                 sigma=matrix(c(1,Inf,1,1,1,1,Inf,1),2,4),
                 tau=0,N=c(2,3),psi=5,gamma=9)
construct_CovMat(SumCl=2,timepoints=4,
                 sigma=matrix(c(1,Inf,1,1,1,1,Inf,1),2,4),
                 tau  =matrix(c(0,1)),
                 N=c(2,3),psi=5,gamma=9)
construct_CovMat(SumCl=2,timepoints=4,
                 sigma=matrix(c(1,Inf,1,1,1,1,Inf,1),2,4),
                 tau  =matrix(c(1,1)),
                 tauAR=c(1,.7),
                 N=c(2,3),psi=5,gamma=9)
construct_CovMat(SumCl=2,timepoints=4,
                 sigma=matrix(c(1,Inf,1,1,1,1,Inf,1),2,4),
                 tau  =matrix(c(1,1)),
                 tauAR=c(1,1),
                 eta    =1,
                 trtMat =matrix(c(0,0,1,0,1,1,1,1),2,4),
                 N=c(2,3),psi=5,gamma=9, INDIV_LVL=TRUE)

construct_CovSubMat(N=2,
                    sigma = rep(sqrt(2)/2,4),
                    tau =rep(1,4),
                    gamma=rep(9,4),
                    psi=rep(5,4),
                    timepoints = 4,
                    INDIV_LVL = FALSE)

construct_CovMat(SumCl=2,timepoints=3,sigma=4,tau=1,N=2,psi=.1,INDIV_LVL = TRUE)
construct_CovMat(SumCl=2,timepoints=3,sigma=4,tau=1,N=2,psi=0,INDIV_LVL = TRUE)

################################################################################
## compute_wlsPower #####

## some tests for denomDF adjustment: #####
compute_wlsPower(DesMat=construct_DesMat(Cl=c(4,2,2,4)),
    EffSize=.5,sigma=1,tau=.1,N=1,dfAdjust="none",sig.level=.05,verbose=T)
compute_wlsPower(DesMat=construct_DesMat(Cl=c(4,2,2,4)),
    EffSize=.5,sigma=1,tau=.1,N=1,dfAdjust="between-within",sig.level=.05,verbose=T)
compute_wlsPower(DesMat=construct_DesMat(Cl=c(4,2,2,4)),
    EffSize=.5,sigma=1,tau=.1,N=1,dfAdjust="containment",sig.level=.05,verbose=T)
compute_wlsPower(DesMat=construct_DesMat(Cl=c(4,4,4,4)),
    EffSize=.5,sigma=1,tau=.1,N=1,dfAdjust="between-within",sig.level=.05,verbose=T)



compute_wlsPower(DesMat=construct_DesMat(Cl=c(337,337),dsntype="parallel",timepoints=1),
    EffSize=.25,sigma=1,tau=1,N=1,dfAdjust="none",sig.level=.05,verbose=F)

Cl <- rep(10,10)
DesMat <- construct_DesMat(Cl=Cl)
DesMat_prl <- construct_DesMat(Cl=c(40,40),dsntype="parallel",timepoints=3)


SteppedPower:::compute_wlsPower(DesMat=DesMat, EffSize=.05, sigma=1,tau=.3,N=1,
                 dfAdjust="none",sig.level=.05,verbose=F)
wlsPower(DesMat=DesMat,mu0=0,mu1=.05,sigma=1,tau=.3,verbose=F)



SteppedPower:::compute_wlsPower(DesMat=DesMat_prl, EffSize=.5, sigma=1,tau=.3,N=1,
                 dfAdjust="none",sig.level=.05,verbose=F)
wlsPower(DesMat=DesMat_prl,mu0=0,mu1=.5,sigma=1,tau=.3,verbose=F)
wlsPower(Cl=c(4,4),timepoints=3,dsntype="parallel",
              mu0=0,mu1=.5,sigma=1,tau=.3,verbose=F)





DesMat <- construct_DesMat(rep(1,8))
wlsPower(DesMat=DesMat,mu0=0,mu1=.05,sigma=1,tau=.3,N=720,verbose=F)
wlsPower(DesMat=DesMat,mu0=0,mu1=.05,sigma=1,tau=.3,Power=.9,verbose=F)

uniroot(function(N){compute_wlsPower(DesMat,EffSize=.05,sigma=1,tau=.3,N=N)$Power-.9},
        interval=c(1,1000))

wlsPower(DesMat=DesMat_prl,mu0=0,mu1=.15,sigma=1,tau=.15,N=17,verbose=F)
wlsPower(DesMat=DesMat_prl,mu0=0,mu1=.15,sigma=1,tau=.15,Power=.9,verbose=F,N_range = c(1,2000))


################################################################################
## wlsPower #####
wlsPower(mu0=0,mu1=.1,sigma=1,tau=.01,Cl=c(1,1,1,1),verbose=F)

Cl <- c(2,2,2,2,2); sigma <- 4; tau <- 0.9 ; mu0=0 ; mu1=1.5
wlsPower(Cl=Cl,sigma=sigma,tau=tau,mu0=mu0,mu1=mu1,verbose=F)
swCRTdesign::swPwr(swCRTdesign::swDsn(Cl),distn="gaussian",1,0,EffSize,
      tau=tau,eta=0,rho=0,sigma=sigma)

wlsPower(mu0=0,mu1=1,Cl=c(1,1,1,1,1),sigma=2 ,        tau=0.2, N=c(1,1,1,1,1) )
wlsPower(mu0=0,mu1=1,Cl=c(1,1,1,1,1),sigma=2*sqrt(2) ,tau=0.2, N=c(2,2,2,2,2) )

Cl=c(3,3) ; mu0=0 ; mu1=1 ; sigma=1 ; tau=0.1 ; family=gaussian() ; timepoints=1 ;
dsntype="parallel" ; N <- NULL ; sig.level=0.05
wlsPower(Cl=c(5,5),mu0=0,mu1=.2,sigma=1,tau=0.1,timepoints=2,dsntype="parallel")

microbenchmark::microbenchmark(
wls_parBl1 <- wlsPower(mu0=0,mu1=0.1,sigma=1,tau=1,Cl=c(50,50),
                           dsntype="parallel_baseline",timepoints=51)
,
wls_parBl2 <- wlsPower(mu0=0,mu1=0.1,sigma=1,tau=1,Cl=c(50,50),
                           dsntype="parallel_baseline",timepoints=c(15,36))
,
wls_swd    <- wlsPower(mu0=0,mu1=0.1,sigma=1,tau=1,Cl=rep(2,50),
                            dsntype="SWD")
,times=10)
system.time(sw <- swCRTdesign::swPwr(swCRTdesign::swDsn(rep(2,50)),distn="gaussian",
                               n=1,mu0=0,mu1=.1,
                               sigma=1,tau=1,eta=0,rho=0,gamma=0))
wls_parBl1$Power
wls_parBl2$Power
wls_swd$Power

Rprof(interval=0.005,line.profiling=TRUE)
wls_swd    <- wlsPower(mu0=0,mu1=0.1,sigma=1,tau=1,Cl=rep(3,100),
                            dsntype="SWD")
Rprof(NULL)
summaryRprof()


wlsPower(mu0=0,mu1=.1,sigma=1,tau=.3,Cl=c(2,2,2,2,2),Power=.9,verbose=F)
wlsPower(mu0=0,mu1=.1,sigma=1,tau=.3,Cl=c(2,2,2,2,2),N=224,verbose=F)
swCRTdesign::swPwr(swCRTdesign::swDsn(c(2,2,2,2,2)),distn="gaussian",n=224,
                   mu0=0,mu1=.1,tau=.3,eta=0,rho=0,gamma=0,sigma=1)

wlsPower(mu0=0,mu1=.1,sigma=1,tau=.3,eta=.1,Cl=c(2,2,2,2,2),N=224)
swCRTdesign::swPwr(swCRTdesign::swDsn(c(2,2,2,2,2)),distn="gaussian",n=224,
                   mu0=0,mu1=.1,tau=.3,eta=.1,rho=0,gamma=0,sigma=1)

cls <- rep(5,10)
microbenchmark::microbenchmark(
a1 <- wlsPower(mu0=0,mu1=.1,sigma=1,tau=.3,eta=.1,rho=1,Cl=cls,N=224)
,
a2 <- swCRTdesign::swPwr(swCRTdesign::swDsn(cls),distn="gaussian",n=1,
                   mu0=0,mu1=.1,tau=0,eta=0,rho=0,gamma=0,sigma=5,retDATA=TRUE)
,times=10)
a2$Wmat
a1$CovarianceMatrix[1:6,1:6];a2$Wmat[1:6,1:6]
a1$Power ; a2$pwrWLS

wlsPower(mu0=0,mu1=.1,sigma=1,tau=.3,eta=.1,rho=1,Cl=c(2,2,2,2,2),N=224,timeAdjust="none")


wlsPower(mu0=0,mu1=.05,sigma=1,tau=.3,Cl=rep(2,20),Power=.9,verbose=F)
wlsPower(mu0=0,mu1=.05,sigma=1,tau=.3,Cl=rep(2,20),N=60,verbose=F)

wlsPower(mu0=0,mu1=.02,sigma=1,tau=.0,
              Cl=c(50,50),timepoints=5,dsntype="parallel",
              Power=.9,verbose=F)
wlsPower(mu0=0,mu1=.02,sigma=1,tau=.0,
              Cl=c(50,50),timepoints=5,dsntype="parallel",
              N=300,verbose=F)


wlsPower(Cl=c(1,1,1),trtDelay=c(.3,.7),timeAdjust="None",mu0=0,mu1=.1,sigma=1,tau=.1,verbose=T)
wlsPower(Cl=c(1,1,1,1),trtDelay=c(.3,.7),timeAdjust="None",mu0=0,mu1=.1,sigma=1,tau=.1,Power=.8,verbose=T)
wlsPower(Cl=c(1,1,1,1),trtDelay=c(.3,.7),timeAdjust="factor",mu0=0,mu1=.1,sigma=1,tau=.01,Power=.8,verbose=T)

wlsPower(Cl=c(2,2,2,2,2,2),trtDelay=c(.3,.7),timeAdjust="None",mu0=0,mu1=.1,sigma=1,tau=.01,Power=.8,verbose=T)
wlsPower(Cl=c(2,2,2,2,2,2),trtDelay=c(.3,.7),timeAdjust="factor",mu0=0,mu1=.1,sigma=1,tau=.01,Power=.8,verbose=T)


wlsPower(Cl=c(2,2,2,0,2,2,2,0),mu0=0,mu1=.01,
              sigma=sqrt(.025*.975), timeAdjust="none",
              tau=0.00254,trtDelay=.5, N=58, verbose=2)
wlsPower(Cl=c(2,2,2,0,2,2,2,0),mu0=0,mu1=.01,
              sigma=sqrt(.025*.975), timeAdjust="periodic", period=4,
              tau=0.00254,trtDelay=.5, N=58, verbose=2)

a <- wlsPower(Cl=c(2,2,2,0,2,2,2),mu0=0,mu1=.01, sigma=sqrt(.025*.975),
              tau=0.0254, gamma=15, N=58, verbose=2)
a$CovarianceMatrix


#### tauAR swd ####
cl <- rep(6,3) ; si <- 2 ; tau <- .1 ; ES <- 1
mod1 <- wlsPower(Cl=cl, sigma=si, tau=tau, tauAR=NULL, mu0=0,mu1=ES, verbose=2)
mod2 <- wlsPower(Cl=cl, sigma=si, tau=tau, tauAR=1,    mu0=0,mu1=ES, verbose=2)
mod3 <- wlsPower(Cl=cl, sigma=si, tau=tau, tauAR=.7,   mu0=0,mu1=ES, verbose=2)
mod4 <- wlsPower(Cl=cl, sigma=si, tau=tau, tauAR=.2,   mu0=0,mu1=ES, verbose=2)
mod5 <- wlsPower(Cl=cl, sigma=si, tau=tau, tauAR=.1,   mu0=0,mu1=ES, verbose=2)
mod6 <- wlsPower(Cl=cl, sigma=si, tau=tau, tauAR=0,    mu0=0,mu1=ES, verbose=2)
mod7 <- wlsPower(Cl=cl, sigma=sqrt(si^2+tau^2),        mu0=0,mu1=ES, verbose=2)
mod8 <- wlsPower(Cl=cl, sigma=si, tau=0,   gamma=tau,  mu0=0,mu1=ES, verbose=2)

mod1$CovarianceMatrix[1:4,1:4]
mod2$CovarianceMatrix[1:4,1:4]
mod3$CovarianceMatrix[1:4,1:4]
mod4$CovarianceMatrix[1:4,1:4]
mod5$CovarianceMatrix[1:4,1:4]
mod6$CovarianceMatrix[1:4,1:4]
mod7$CovarianceMatrix[1:4,1:4]
mod8$CovarianceMatrix[1:4,1:4]

mod1$Power
mod2$Power
mod3$Power
mod4$Power
mod5$Power
mod6$Power
mod7$Power
mod8$Power

#### tauAR parallel ####
cl <- c(10,10) ; si <- 2 ; tau <- 2 ; ES <- 1 ; tp <- 6

#debugonce(wlsPower)
#debugonce(compute_wlsPower)

mod1 <- wlsPower(Cl=cl, timepoints=tp, sigma=si, tau=tau, tauAR=NULL, mu0=0,mu1=ES, dsntype="parallel", verbose=2)
mod2 <- wlsPower(Cl=cl, timepoints=tp, sigma=si, tau=tau, tauAR=1,    mu0=0,mu1=ES, dsntype="parallel", verbose=2)
mod3 <- wlsPower(Cl=cl, timepoints=tp, sigma=si, tau=tau, tauAR=.8,   mu0=0,mu1=ES, dsntype="parallel", verbose=2)
mod4 <- wlsPower(Cl=cl, timepoints=tp, sigma=si, tau=tau, tauAR=.5,   mu0=0,mu1=ES, dsntype="parallel", verbose=2)
mod5 <- wlsPower(Cl=cl, timepoints=tp, sigma=si, tau=tau, tauAR=0,    mu0=0,mu1=ES, dsntype="parallel", verbose=2)
mod6 <- wlsPower(Cl=cl, timepoints=tp, sigma=sqrt(si^2+tau^2),        mu0=0,mu1=ES, dsntype="parallel", verbose=2)
mod7 <- wlsPower(Cl=cl, timepoints=tp, sigma=si, tau=0,   gamma=tau,  mu0=0,mu1=ES, dsntype="parallel", verbose=2)


mod1$CovarianceMatrix[1:tp,1:tp]
mod2$CovarianceMatrix[1:tp,1:tp]
mod3$CovarianceMatrix[1:tp,1:tp]
mod4$CovarianceMatrix[1:tp,1:tp]
mod5$CovarianceMatrix[1:tp,1:tp]
mod6$CovarianceMatrix[1:tp,1:tp]
mod7$CovarianceMatrix[1:tp,1:tp]

mod1$Power
mod2$Power
mod3$Power
mod4$Power
mod5$Power
mod6$Power


#### incomplete Stepped Wedge design  ####
Des_SWD22 <- construct_DesMat(Cl=c(2,2))

#debugonce(wlsPower)
a <- wlsPower(Cl=c(2,2,2), mu=0, mu1=1, sigma=1, tau=2, incomplete=2,verbose = T)
plot(a)

IncMat1 <- matrix(c(1,1,0,1,1,1,1,1,1,0,1,1),3,4)
b <- wlsPower(Cl=c(2,2,2), mu=0, mu1=1, sigma=1, tau=2,
                   incomplete=IncMat1,verbose = T)
plot(b)
all.equal(a,b)

IncMat2 <- IncMat1[rep(1:3,each=2),]

c <- wlsPower(Cl=c(2,2,2), mu=0, mu1=1, sigma=1, tau=2,
                   incomplete=IncMat2,verbose = T)
plot(c)
all.equal(a,c)

## closed cohort ####

# debugonce(wlsPower)
# debugonce(compute_wlsPower)
# debugonce(construct_CovMat)

wlsPower(Cl=c(2,3),mu=0,mu1=1,tau=.1)
wlsPower(Cl=c(2,3),mu=0,mu1=1,tau=.1,psi=0,N=1)

sig <- .05
wlsPower(Cl=c(2,3),mu=0,mu1=.05,sigma=sig, tau=10,N=3:7)
wlsPower(Cl=c(2,3),mu=0,mu1=.05,sigma=sig, tau=10,psi=0,N=3:7)
wlsPower(Cl=c(2,3),mu=0,mu1=.05,sigma=sig, tau=10,psi=0,N=3:7,INDIV_LVL=TRUE)

wlsPower(Cl=c(2,3),mu=0,mu1=.05,sigma=sig, tau=10,N=matrix(11:25,5))

a  <- wlsPower(Cl=c(4,4),mu=0,mu1=1,tau=.1,verbose = TRUE)
b1 <- wlsPower(Cl=c(1,1),mu=0,mu1=1,tau=0,psi=.1,N=c(4,4), verbose = TRUE)
b2 <- wlsPower(Cl=c(1,1),mu=0,mu1=1,tau=0,psi=.1,N=c(4,4), verbose = TRUE, INDIV_LVL = TRUE)
c1 <- wlsPower(Cl=c(2,2),mu=0,mu1=1,tau=0,psi=.1,N=c(2,2,2,2), verbose=2)
c2 <- wlsPower(Cl=c(2,2),mu=0,mu1=1,tau=0,psi=.1,N=c(2,2,2,2), verbose=2, INDIV_LVL = TRUE)
a;b1;b2;c1;c2

plot(a)
plot(b2)
plot(c2)

d1 <- wlsPower(Cl=c(2,2),mu=0,mu1=1,sigma=.001,tau=1,psi=.1,N=c(2,2,2,2),eta=1, verbose=2)
d2 <- wlsPower(Cl=c(2,2),mu=0,mu1=1,sigma=.001,tau=1,psi=.1,N=c(2,2,2,2),eta=1, verbose=2, INDIV_LVL = TRUE)
d1;d2

e <- wlsPower(Cl=c(2,2,2),mu=0,mu1=1,sigma=1,tau=1,eta=2,psi=.1,N=c(2,2,2,2,2,2), incomplete=2, verbose=2)
plot(e)

f <- wlsPower(Cl=c(2,2,2,2),mu=0,mu1=1,sigma=1,tau=1,eta=.5,gamma=.01, psi=1,N=rep(1:4,2), verbose=2)
f;plot(f)

g <- wlsPower(Cl=rep(2,6),mu=0,mu1=.5,sigma=1,tau=1,eta=.5,gamma=.01, psi=0, N=100, verbose=2)
g;plot(g)
wlsPower(Cl=rep(2,6),mu=0,mu1=.5,sigma=1,tau=1,eta=.5,gamma=.01, N=100, verbose=2)


Cl <- rep(20,20) ; tp <- length(Cl)+1
nInd <- 500
tau <- .2 ; sigma <- .3 ; psi <- .3 ; eta <- 0 ; gamma <- .1

system.time({
  Rprof(filename="Rprof.out", interval=.005, line.profiling=TRUE)
  b <- wlsPower(Cl=Cl, sigma=sigma,
              tau=tau, eta=eta, psi=psi, gamma=gamma,
              mu0=0,mu1=.05, N=nInd)
  Rprof(NULL)
})
summaryRprof(lines="show")

## other design types for closed cohorts ####

nCl <- 10 ; N <- 20
microbenchmark::microbenchmark({
wlsPower(Cl=c(nCl,nCl), timepoints=3, dsntype="parallel", mu0=0, mu1=1,
              sigma=1, tau=2, N=N)
# },{
wlsPower(Cl=c(nCl,nCl), timepoints=3, dsntype="parallel", mu0=0, mu1=1,
              sigma=1, tau=2, psi=0, N=N)
# },{
wlsPower(Cl=c(nCl,nCl), timepoints=3, dsntype="parallel", mu0=0, mu1=1,
              sigma=1, tau=2, psi=0, N=N, INDIV_LVL = TRUE)
# })

wlsPower(Cl=c(nCl,nCl), timepoints=3, dsntype="parallel", mu0=0, mu1=1,
              sigma=1, tau=2, psi=1, N=N)
wlsPower(Cl=c(nCl,nCl), timepoints=3, dsntype="parallel", mu0=0, mu1=1,
              sigma=1, tau=2, psi=1, N=N)


wlsPower(Cl=c(nCl,nCl), timepoints=4, dsntype="parallel_baseline", mu0=0, mu1=1,
              sigma=1, tau=2, psi=1, N=1)
wlsPower(Cl=c(nCl,nCl), timepoints=4, dsntype="parallel_baseline", mu0=0, mu1=1,
              sigma=1, tau=sqrt(5), psi=0, N=1)


wlsPower(Cl=c(nCl,nCl), timepoints=4, dsntype="parallel_baseline", mu0=0, mu1=1,
              sigma=1, tau=2, psi=1, eta=.5, N=1)
wlsPower(Cl=c(nCl,nCl), timepoints=4, dsntype="parallel_baseline", mu0=0, mu1=1,
              sigma=1, tau=2, psi=1, eta=.5, N=1, INDIV_LVL = TRUE)




microbenchmark::microbenchmark( cat(" ") ,  cat(".") , cat("_") )



## ICC and CAC ####

wlsPower(Cl=c(2,2,2), mu0=0, mu1=1, icc=.1, cac=.05)
icc_to_RandEff(.1,.05,1)
wlsPower(Cl=c(2,2,2), mu0=0, mu1=1, tau=.0745, gamma=.325)


################################################################################
## plot.wlsPower #####
plot(wlsPower(Cl=c(2,5),sigma=1,tau=0.1, mu0=0, mu1=1,
                            timepoints=5,dsntype="parallel",verbose=2))


DesMat11 <- construct_DesMat(Cl=c(2,1))
CovMat11 <- construct_CovMat(3,3,1,2)
IncMat11 <- matrix(c(0,1,1,1,1,1,1,1,1),3,3)
IncMat22 <- matrix(c(1,1,1,1,0,1,1,1,1),3,3)
IncMat31 <- matrix(c(1,1,0,1,1,1,1,1,1),3,3)
IncMat2122 <- matrix(c(1,0,1,1,0,1,1,1,1),3,3)

wmp <- wlsPower(DesMat=DesMat11, mu0=0, mu1=1, sigma=1, tau=2,
                          incomplete = NULL, verbose=2)
wmpinc11 <- wlsPower(DesMat=DesMat11, mu0=0, mu1=1, sigma=1, tau=2,
                          incomplete = IncMat11, verbose=2)
wmpinc22 <- wlsPower(DesMat=DesMat11, mu0=0, mu1=1, sigma=1, tau=2,
                          incomplete = IncMat22, verbose=2)
wmpinc31 <- wlsPower(DesMat=DesMat11, mu0=0, mu1=1, sigma=1, tau=2,
                          incomplete = IncMat31, verbose=2)
wmpinc2122 <- wlsPower(DesMat=DesMat11, mu0=0, mu1=1, sigma=1, tau=2,
                          incomplete = IncMat2122, verbose=2)
plot(wmp)
plot(wmpinc11)
plot(wmpinc22)
plot(wmpinc31)

wmp_to_Var <- function(wmp){
  dsnmatrix <- wmp$DesignMatrix$dsnmatrix
  CovMat    <- wmp$CovarianceMatrix

  tmpmat  <- t(dsnmatrix) %*% Matrix::chol2inv(Matrix::chol(CovMat))
  VarMat  <- Matrix::solve(tmpmat %*% dsnmatrix)
  return(VarMat)
}

sqrt(wmp_to_Var(wmp)[1,1])
sqrt(wmp_to_Var(wmpinc11)[1,1])
sqrt(wmp_to_Var(wmpinc22)[1,1])
sqrt(wmp_to_Var(wmpinc31)[1,1])
sqrt(wmp_to_Var(wmpinc2122)[1,1])



dsnmatrix <- wmp$DesignMatrix$dsnmatrix
CovMat    <- wmp$CovarianceMatrix

tmpmat  <- t(dsnmatrix) %*% Matrix::chol2inv(Matrix::chol(CovMat))
VarMat  <- Matrix::solve(tmpmat %*% dsnmatrix)
ProjMat <- matrix((VarMat %*% tmpmat)[1,], nrow=3, byrow=TRUE)

HatMat <- dsnmatrix %*% VarMat %*% tmpmat
hi <- matrix(diag(HatMat),nrow=3,byrow=TRUE)
ProjMat/(1-hi)

## E[DFBETA] ####
i <- 1
V_inv <- solve(CovMat)
v_ii  <- diag(V_inv)[i]
v_i   <- V_inv[,i]
sd_e_i<- sigma
E_h_i <- (1/v_ii)*v_i*sigma
q_i   <-

g_i <-

## Leverage ####

V_inv <- solve(CovMat)
X     <- wmp$DesignMatrix$dsnmatrix
P  <- V_inv %*% X %*% solve(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv
P2 <- V_inv %*% HatMat %*% V_inv
round(P-P2,5)

tmpmat2 <- t(X) %*% solve(CovMat)
round(tmpmat-tmpmat2,5)

HatMat2 <- X %*% solve(tmpmat2 %*% X) %*% tmpmat2
round(HatMat-HatMat2,5)


Q     <- V_inv - P
round(matrix(diag(Q),3,3,byrow=TRUE),2)

############################### plot random effect robustness ##################

library(plotly)

x <- wlsPower(Cl=c(2,2,2,2,1,1,2,2), trtDelay = .5, mu0=0, mu1=.01, N=500,
                   sigma=sqrt(0.0244), tau=0.005, eta=0.005, verbose=2)
x <- wlsPower(Cl=c(8,8), trtDelay = .5, mu0=0, mu1=.01, N=500,
                   sigma=sqrt(0.0244), tau=0.02,
                   eta=0.02, verbose=2)
x$Params


p <- x$Params

Range <- seq(0,2,length.out=4)
rhoRange <- c(-.5,0,.5)
if(!is.null(p$rho)) rhoRange <- append(rhoRange,p$rho)

tauRange <- lapply(Range,function(x) x * p$tau)
etaRange <- lapply(Range,function(x) x * p$eta)

taugrid <- expand.grid(rhoRange=rhoRange,tauRange=tauRange)

tauPwr   <-  mapply(wlsPower,
                    tau=tauRange,
                    # tau=as.list(taugrid$tauRange),
                    # rho=as.list(taugrid$rhoRange),
                    MoreArgs = list(DesMat=x$DesignMatrix,
                                    N=p$N,
                                    sigma=p$sigma,
                                    eta=p$eta,
                                    mu0=p$mu0,
                                    mu1=p$mu1,
                                    verbose=0),
                    SIMPLIFY = TRUE)

# plot_ly(type="scatter",mode="lines+marker",
#         x=~tauRange, y=~tauPwr, split=~rhoRange,
#         data=cbind(taugrid,tauPwr))

etaPwr   <-  mapply(wlsPower,eta=as.list(etaRange),
                    MoreArgs = list(DesMat=x$DesignMatrix,
                                    N=p$N,
                                    sigma=p$sigma,
                                    tau=p$tau,
                                    mu0=p$mu0,
                                    mu1=p$mu1,
                                    verbose=0),
                    SIMPLIFY = TRUE)


tauplt <- plotly::plot_ly()
for(i in 1:1){
  tauplt <- add_lines(tauplt,x=Range, y=tauPwr)
}

etaplt <- plotly::plot_ly()
for(i in 1:1){
  etaplt <- add_lines(etaplt,x=Range, y=etaPwr)
}



add_original_Pwr <- function(plt,xRange,pwr){
  plt <- add_segments(plt,
                      x=xRange[1], xend=xRange[2],
                      y=pwr, yend=pwr,
                      line=list(color='darkgrey', dash="dash", width=5))
}



plot_randeffs <- function(X,Range=seq(0,2,length.out=10)){

  p        <- X$Params
  tauRange <- lapply(Range, function(x) x * p$tau)
  etaRange <- lapply(Range, function(x) x * p$eta)
  tauPwr   <-  mapply(wlsPower,
                      tau=as.list(tauRange),
                      MoreArgs = list(DesMat=X$DesignMatrix,
                                      N=p$N,
                                      sigma=p$sigma,
                                      eta=p$eta,
                                      mu0=p$mu0,
                                      mu1=p$mu1,
                                      verbose=0) )

  etaPwr   <-  mapply(wlsPower,
                      eta=as.list(etaRange),
                      MoreArgs = list(DesMat=X$DesignMatrix,
                                      N=p$N,
                                      sigma=p$sigma,
                                      tau=p$tau,
                                      mu0=p$mu0,
                                      mu1=p$mu1,
                                      verbose=0) )

  tauplt <- plotly::plot_ly()
  tauplt <- add_lines(tauplt,x=Range, y=tauPwr)
  tauplt <- add_original_Pwr(tauplt,xRange=range(Range),pwr=x$Power)

  etaplt <- plotly::plot_ly()
  etaplt <- add_lines(etaplt,x=Range, y=etaPwr)
  etaplt <- add_original_Pwr(etaplt,xRange=range(Range),pwr=x$Power)

  subplot(  tauplt, etaplt, shareY=TRUE  )
}


x <- wlsPower(Cl=c(2,2,2,2), trtDelay=.3, sigma=1, tau=.2, eta=.3,
                   mu0=0, mu1=1, verbose = 2)

plot_randeffs(x)


plot_wls3(x)

plot_wls3 <- function(x){
  p      <- x$Params
  tauRange <- seq(0, 2*p$tau,length.out=5)
  etaRange <- if(is.null(p$eta)) seq(0,max(tauRange),length.out=5) else seq(0,2*p$eta,length.out=5)

  grid <- expand.grid(tau=tauRange,eta=etaRange)
  pwr <- list()
  for(i in 1:dim(grid)[1]){
    pwr[[i]] <- wlsPower(DesMat=x$DesignMatrix,mu0=0,mu1=.01,sigma=p$sigma,N=200,
                              tau=grid$tau[i], eta=grid$eta[i])[["Power"]]
  }
  tmp <- data.frame(grid,pwr=unlist(pwr))
  plot_ly(type="scatter",mode="lines+marker",x=~tau, y=~pwr, split=~eta, data=tmp)
}


################################################################################

## compare_designs #####
compare_designs(mu0=0,mu1=1, sigma=1 ,tau=.3, Cl=c(2,2,2,2))
compare_designs(mu0=0,mu1=1, sigma=1 ,tau=.5, Cl=c(2,2,2,2))
compare_designs(mu0=0,mu1=1, sigma=1 ,tau=.7, Cl=c(2,2,2,2))
compare_designs(mu0=0,mu1=1, sigma=1 ,tau=1 , Cl=c(2,2,2,2))


## tTestPwr scaled Wald-test #####

# tTestPwr <- function(d,se,df,sig.level=0.05){
  dsz <- abs(d/se)
  q   <- qt(sig.level/2, df=df, lower=FALSE)
  Pwr <- pt(dsz + q, df=df) + pt(-dsz + q, df=df)
#  return(Pwr)
# }

d <- .6; se <- 1; df <- 18
# debugonce(tTestPwr)
# debugonce(tTestPwr2)
# debugonce(pwr.t.test)
SteppedPower::tTestPwr(.6,1,Inf)
pwr::pwr.t.test(n=10,d=.6,sig.level=.05)
tTestPwr2(.6,1,18)
tTestPwr2(.6,1,Inf)
View(pwr::pwr.t.test)
View(tTestPwr)




################################################################################

library(swCRTdesign)
View(swPwr)
dsn <- swDsn(c(2,2,0,2))
dsn
construct_DesMat(c(2,2,0,2))
## convert swDsn to DesMat


################################################################################
## alpha0_1_2 to random effects

# alpha0 := within-period corr of different individuals
# alpha1 := different indivs across periods
# alpha2 := same individual across periods

sigSq <- 1
tauSq <- 2

alpha_0 <- 1/2
alpha_1 <- 1/3

tauSq   <- (sigSq*alpha_1)/(1-alpha_1)
gammaSq <-

################################################################################

debugonce(muCond_to_muMarg)
muMarg_to_muCond(.4,1)
muMarg_to_muCond(.02,.42)
muMarg_to_muCond(.0001,.42)

muCond_to_muMarg(.001,.003)
muCond_to_muMarg(.01,0.002)

## PLOTS for muMarg_to_muCond ####


mus  <- seq(1e-5,1-1e-5, 5e-4)
taus <- seq(1e-5,2 , 5e-4)
y1 <- sapply(mus, muMarg_to_muCond, tauLin=.4)  ## changing mu, fixed    tau
y2 <- sapply(taus,muMarg_to_muCond, muMarg=.001)## fixed    mu, changing tau
y3 <- mapply(muMarg_to_muCond, muMarg=mus, tauLin=3*sqrt(mus*(1-mus)))
y4 <- mapply(muMarg_to_muCond, muMarg=mus, tauLin=4*sqrt(mus*(1-mus)))
y5 <- mapply(muMarg_to_muCond, muMarg=mus, tauLin=5*sqrt(mus*(1-mus)))

plot(mus,y1, "l")
lines(mus,mus,col=2)

plot(taus,y2, "l")

plot(mus, (y3-mus),"l")
plot(mus,y3, "l")
lines(mus,mus,col=2)
lines(mus,y4,col=3)
lines(mus,y5,col=4)
maxDiff <- which.max(y3-mus)
mus[maxDiff]  ; y3[maxDiff]


## PLOTS for muCond_to_muMarg ####

y1 <- sapply(mus, muCond_to_muMarg, tauLin=.4)
y2 <- sapply(taus,muCond_to_muMarg, muCond=.25)

plot(mus,y1, type="l")
lines(mus,mus,col=2)
plot(mus,(y1-mus),"l")

plot(taus,y2/.25, "l")
