
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
construct_trtMat(Cl=c(1,1,1), trtDelay=NULL, design="SWD",      timepoints=4)
construct_trtMat(Cl=c(1,1),   trtDelay=NULL, design="parallel", timepoints=4)
construct_trtMat(Cl=c(1,1,1), trtDelay=c(.2,.4),  design="SWD",      timepoints=4)
construct_trtMat(Cl=c(1,1),   trtDelay=c(.2,.4),  design="parallel", timepoints=4)
construct_trtMat(Cl=c(2,2),   trtDelay=NULL,  design="parallel_baseline", timepoints=3)
construct_trtMat(Cl=c(2,2),   trtDelay=NULL,  design="crossover", timepoints=c(2,2))
# debugonce(construct_trtMat)

## costruct_DesMat #####
construct_DesMat(Cl=c(2,0,1))
Des_PB <- construct_DesMat(Cl=c(2,2))

construct_DesMat(Cl=c(5,5),design="parallel",timepoints=NULL)
construct_DesMat(Cl=c(5,5),design="parallel",timepoints=1)
Des_PB <-construct_DesMat(Cl=c(2,2),design="parallel",timepoints=3)

construct_DesMat(Cl=c(2,2),design="parallel_baseline",timepoints=NULL)
construct_DesMat(Cl=c(2,2),design="parallel_baseline",timepoints=2)
Des_PB <- construct_DesMat(Cl=c(2,2),design="parallel_baseline",timepoints=c(1,2))

construct_DesMat(Cl=c(1,1,1),trtDelay=c(.3,.7))
construct_DesMat(Cl=c(1,1),  trtDelay=c(.3,.7),design="parallel")
construct_DesMat(Cl=c(1,1),  trtDelay=c(.3,.7),design="parallel_baseline")

construct_DesMat(Cl=c(1,1,1),trtDelay=c(.3,.7),timeAdjust="none")
construct_DesMat(Cl=c(1,1),trtDelay=c(.3,.7),design="parallel_baseline",timeAdjust="linear")

construct_DesMat(Cl=c(1,1,1,1),timeAdjust="factor")
construct_DesMat(Cl=c(1,1,1,1),timeAdjust="periodic")
construct_DesMat(Cl=rep(1,6),timeAdjust="periodic",period=4)

construct_DesMat(Cl=c(1,1,2),trtDelay=.5,timeBlk=diag(4))

################################################################################
## CovBlk #####
timepoints=3; sigma=2; tau=rep(2,3)
construct_CovBlk(timepoints=5, sigma=2, tau=2)
construct_CovBlk(timepoints=3, sigma=2, tau=rep(2,3))
construct_CovBlk(timepoints=3, sigma=c(1,2,3), tau=2)

eta <- c(0,.2,.2)
construct_CovBlk(timepoints,sigma,tau,eta)

rho <- .3
construct_CovBlk(timepoints,sigma,tau,eta,rho)

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
construct_CovMat(SumCl=3,timepoints=4,sigma=sqrt(10),tau=1,eta=sqrt(.1),rho=.5,trtMat=trtMat)


construct_CovMat(SumCl=3, timepoints=4, sigma=3, tau=2, tauAR=.8)

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



compute_wlsPower(DesMat=construct_DesMat(Cl=c(337,337),design="parallel",timepoints=1),
    EffSize=.25,sigma=1,tau=1,N=1,dfAdjust="none",sig.level=.05,verbose=F)

Cl <- rep(10,10)
DesMat <- construct_DesMat(Cl=Cl)
DesMat_prl <- construct_DesMat(Cl=c(40,40),design="parallel",timepoints=3)


SteppedPower:::compute_wlsPower(DesMat=DesMat, EffSize=.05, sigma=1,tau=.3,N=1,
                 dfAdjust="none",sig.level=.05,verbose=F)
wlsMixedPower(DesMat=DesMat,mu0=0,mu1=.05,sigma=1,tau=.3,verbose=F)



SteppedPower:::compute_wlsPower(DesMat=DesMat_prl, EffSize=.5, sigma=1,tau=.3,N=1,
                 dfAdjust="none",sig.level=.05,verbose=F)
wlsMixedPower(DesMat=DesMat_prl,mu0=0,mu1=.5,sigma=1,tau=.3,verbose=F)
wlsMixedPower(Cl=c(4,4),timepoints=3,design="parallel",
              mu0=0,mu1=.5,sigma=1,tau=.3,verbose=F)





DesMat <- construct_DesMat(rep(1,8))
wlsMixedPower(DesMat=DesMat,mu0=0,mu1=.05,sigma=1,tau=.3,N=720,verbose=F)
wlsMixedPower(DesMat=DesMat,mu0=0,mu1=.05,sigma=1,tau=.3,Power=.9,verbose=F)

uniroot(a <- function(N){compute_wlsPower(DesMat,EffSize=.05,sigma=1,tau=.3,N=N)$Power-.9},
        interval=c(1,1000))

wlsMixedPower(DesMat=DesMat_prl,mu0=0,mu1=.15,sigma=1,tau=.15,N=14,verbose=F)
wlsMixedPower(DesMat=DesMat_prl,mu0=0,mu1=.15,sigma=1,tau=.15,Power=.9,verbose=F,N_range = c(1,20))


################################################################################
## wlsMixedPower #####
wlsMixedPower(mu0=0,mu1=.1,sigma=1,tau=.01,Cl=c(1,1,1,1),verbose=F)

Cl <- c(2,2,2,2,2); sigma <- 4; tau <- 0.9 ; mu0=0 ; mu1=1.5
wlsMixedPower(Cl=Cl,sigma=sigma,tau=tau,mu0=mu0,mu1=mu1,verbose=F)
swCRTdesign::swPwr(swCRTdesign::swDsn(Cl),distn="gaussian",1,0,EffSize,
      tau=tau,eta=0,rho=0,sigma=sigma)

wlsMixedPower(mu0=0,mu1=1,Cl=c(1,1,1,1,1),sigma=2 ,        tau=0.2, N=c(1,1,1,1,1) )
wlsMixedPower(mu0=0,mu1=1,Cl=c(1,1,1,1,1),sigma=2*sqrt(2) ,tau=0.2, N=c(2,2,2,2,2) )

Cl=c(3,3) ; mu0=0 ; mu1=1 ; sigma=1 ; tau=0.1 ; family=gaussian() ; timepoints=1 ;
design="parallel" ; N <- NULL ; sig.level=0.05
wlsMixedPower(Cl=c(5,5),mu0=0,mu1=.2,sigma=1,tau=0.1,timepoints=2,design="parallel")

microbenchmark::microbenchmark(
wls_parBl1 <- wlsMixedPower(mu0=0,mu1=0.1,sigma=1,tau=1,Cl=c(50,50),
                           design="parallel_baseline",timepoints=51)
,
wls_parBl2 <- wlsMixedPower(mu0=0,mu1=0.1,sigma=1,tau=1,Cl=c(50,50),
                           design="parallel_baseline",timepoints=c(15,36))
,
wls_swd    <- wlsMixedPower(mu0=0,mu1=0.1,sigma=1,tau=1,Cl=rep(2,50),
                            design="SWD")
,times=10)
system.time(sw <- swCRTdesign::swPwr(swCRTdesign::swDsn(rep(2,50)),distn="gaussian",
                               n=1,mu0=0,mu1=.1,
                               sigma=1,tau=1,eta=0,rho=0,gamma=0))
wls_parBl1$Power
wls_parBl2$Power
wls_swd$Power


wlsMixedPower(mu0=0,mu1=.1,sigma=1,tau=.3,Cl=c(2,2,2,2,2),Power=.9,verbose=F)
wlsMixedPower(mu0=0,mu1=.1,sigma=1,tau=.3,Cl=c(2,2,2,2,2),N=224,verbose=F)
swCRTdesign::swPwr(swCRTdesign::swDsn(c(2,2,2,2,2)),distn="gaussian",n=224,
                   mu0=0,mu1=.1,tau=.3,eta=0,rho=0,gamma=0,sigma=1)

wlsMixedPower(mu0=0,mu1=.1,sigma=1,tau=.3,eta=.1,Cl=c(2,2,2,2,2),N=224)
swCRTdesign::swPwr(swCRTdesign::swDsn(c(2,2,2,2,2)),distn="gaussian",n=224,
                   mu0=0,mu1=.1,tau=.3,eta=.1,rho=0,gamma=0,sigma=1)

cls <- rep(5,10)
microbenchmark::microbenchmark(
a1 <- wlsMixedPower(mu0=0,mu1=.1,sigma=1,tau=.3,eta=.1,rho=1,Cl=cls,N=224)
,
a2 <- swCRTdesign::swPwr(swCRTdesign::swDsn(cls),distn="gaussian",n=1,
                   mu0=0,mu1=.1,tau=0,eta=0,rho=0,gamma=0,sigma=5,retDATA=TRUE)
,times=10)
a2$Wmat
a1$CovarianceMatrix[1:6,1:6];a2$Wmat[1:6,1:6]
a1$Power ; a2$pwrWLS

wlsMixedPower(mu0=0,mu1=.1,sigma=1,tau=.3,eta=.1,rho=1,Cl=c(2,2,2,2,2),N=224,timeAdjust="none")


wlsMixedPower(mu0=0,mu1=.05,sigma=1,tau=.3,Cl=rep(2,20),Power=.9,verbose=F)
wlsMixedPower(mu0=0,mu1=.05,sigma=1,tau=.3,Cl=rep(2,20),N=60,verbose=F)

wlsMixedPower(mu0=0,mu1=.02,sigma=1,tau=.0,
              Cl=c(50,50),timepoints=5,design="parallel",
              Power=.9,verbose=F)
wlsMixedPower(mu0=0,mu1=.02,sigma=1,tau=.0,
              Cl=c(10,10),timepoints=5,design="parallel",
              N=1000,verbose=F)


wlsMixedPower(Cl=c(1,1,1),trtDelay=c(.3,.7),timeAdjust="None",mu0=0,mu1=.1,sigma=1,tau=.1,verbose=T)
wlsMixedPower(Cl=c(1,1,1,1),trtDelay=c(.3,.7),timeAdjust="None",mu0=0,mu1=.1,sigma=1,tau=.1,Power=.8,verbose=T)
wlsMixedPower(Cl=c(1,1,1,1),trtDelay=c(.3,.7),timeAdjust="factor",mu0=0,mu1=.1,sigma=1,tau=.01,Power=.8,verbose=T)

wlsMixedPower(Cl=c(2,2,2,2,2,2),trtDelay=c(.3,.7),timeAdjust="None",mu0=0,mu1=.1,sigma=1,tau=.01,Power=.8,verbose=T)
wlsMixedPower(Cl=c(2,2,2,2,2,2),trtDelay=c(.3,.7),timeAdjust="factor",mu0=0,mu1=.1,sigma=1,tau=.01,Power=.8,verbose=T)


wlsMixedPower(Cl=c(2,2,2,0,2,2,2,0),mu0=0,mu1=.01,
              sigma=sqrt(.025*.975), timeAdjust="none",
              tau=0.00254,trtDelay=.5, N=58, verbose=TRUE)
wlsMixedPower(Cl=c(2,2,2,0,2,2,2,0),mu0=0,mu1=.01,
              sigma=sqrt(.025*.975), timeAdjust="periodic", period=4,
              tau=0.00254,trtDelay=.5, N=58, verbose=TRUE)

a <- wlsMixedPower(Cl=c(2,2,2,0,2,2,2),mu0=0,mu1=.01, sigma=sqrt(.025*.975),
              tau=0.0254, gamma=15, N=58, verbose=TRUE)
a$CovarianceMatrix


#### gamma #####

CM  <- construct_CovMat(SumCl=2, timepoints=3, sigma=3, tau=sqrt(.1))
CMg <- construct_CovMat(SumCl=2, timepoints=3, sigma=3, tau=sqrt(.1), gamma=100)
DM  <- construct_DesMat(rep(1,2), timeAdjust="none")

Matrix::solve(CM) ; Matrix::solve(CMg)
t(DM$dsnmatrix) %*% Matrix::solve(CM) ; t(DM$dsnmatrix) %*% Matrix::solve(CMg)

tmpmat <- t(DM$dsnmatrix) %*% Matrix::solve(CM)
VarMat <- Matrix::solve(tmpmat %*% DM$dsnmatrix)
VarMat[1,1]

tmpmat <- t(DM$dsnmatrix) %*% Matrix::solve(CMg)
VarMat <- Matrix::solve(tmpmat %*% DM$dsnmatrix)
VarMat[1,1]


#### tauAR swd ####
cl <- rep(6,3) ; si <- 2 ; tau <- .1 ; ES <- 1
mod1 <- wlsMixedPower(Cl=cl, sigma=si, tau=tau, tauAR=NULL, mu0=0,mu1=ES, verbose=TRUE)
mod2 <- wlsMixedPower(Cl=cl, sigma=si, tau=tau, tauAR=1,    mu0=0,mu1=ES, verbose=TRUE)
mod3 <- wlsMixedPower(Cl=cl, sigma=si, tau=tau, tauAR=.7,   mu0=0,mu1=ES, verbose=TRUE)
mod4 <- wlsMixedPower(Cl=cl, sigma=si, tau=tau, tauAR=.2,   mu0=0,mu1=ES, verbose=TRUE)
mod5 <- wlsMixedPower(Cl=cl, sigma=si, tau=tau, tauAR=.1,   mu0=0,mu1=ES, verbose=TRUE)
mod6 <- wlsMixedPower(Cl=cl, sigma=si, tau=tau, tauAR=0,    mu0=0,mu1=ES, verbose=TRUE)
mod7 <- wlsMixedPower(Cl=cl, sigma=sqrt(si^2+tau^2),        mu0=0,mu1=ES, verbose=TRUE)

mod1$CovarianceMatrix[1:4,1:4]
mod2$CovarianceMatrix[1:4,1:4]
mod3$CovarianceMatrix[1:4,1:4]
mod4$CovarianceMatrix[1:4,1:4]
mod5$CovarianceMatrix[1:4,1:4]
mod6$CovarianceMatrix[1:4,1:4]
mod7$CovarianceMatrix[1:4,1:4]

mod1$Power
mod2$Power
mod3$Power
mod4$Power
mod5$Power
mod6$Power
mod7$Power

#### tauAR parallel ####
cl <- c(10,10) ; si <- 2 ; tau <- 2 ; ES <- 1 ; tp <- 6

mod1 <- wlsMixedPower(Cl=cl, timepoints=tp, sigma=si, tau=tau, tauAR=NULL, mu0=0,mu1=ES, design="parallel", verbose=TRUE)
mod2 <- wlsMixedPower(Cl=cl, timepoints=tp, sigma=si, tau=tau, tauAR=1,    mu0=0,mu1=ES, design="parallel", verbose=TRUE)
mod3 <- wlsMixedPower(Cl=cl, timepoints=tp, sigma=si, tau=tau, tauAR=.8,   mu0=0,mu1=ES, design="parallel", verbose=TRUE)
mod4 <- wlsMixedPower(Cl=cl, timepoints=tp, sigma=si, tau=tau, tauAR=.5,   mu0=0,mu1=ES, design="parallel", verbose=TRUE)
mod5 <- wlsMixedPower(Cl=cl, timepoints=tp, sigma=si, tau=tau, tauAR=0,    mu0=0,mu1=ES, design="parallel", verbose=TRUE)
mod6 <- wlsMixedPower(Cl=cl, timepoints=tp, sigma=sqrt(si^2+tau^2),        mu0=0,mu1=ES, design="parallel", verbose=TRUE)

mod1$CovarianceMatrix[1:tp,1:tp]
mod2$CovarianceMatrix[1:tp,1:tp]
mod3$CovarianceMatrix[1:tp,1:tp]
mod4$CovarianceMatrix[1:tp,1:tp]
mod5$CovarianceMatrix[1:tp,1:tp]
mod6$CovarianceMatrix[1:tp,1:tp]

mod1$Power
mod2$Power
mod3$Power
mod4$Power
mod5$Power
mod6$Power


#### incomplete Stepped Wedge design  ####
Des_SWD22 <- construct_DesMat(Cl=c(2,2))

debugonce(wlsMixedPower)
a <- wlsMixedPower(Cl=c(2,2,2), mu=0, mu1=1, sigma=1, tau=2, incomplete=2,verbose = T)
plot(a)

IncMat1 <- matrix(c(1,1,0,1,1,1,1,1,1,0,1,1),3,4)
b <- wlsMixedPower(Cl=c(2,2,2), mu=0, mu1=1, sigma=1, tau=2,
                   incomplete=IncMat1,verbose = T)
plot(b)
all.equal(a,b)

IncMat2 <- IncMat1[rep(1:3,each=2),]

c <- wlsMixedPower(Cl=c(2,2,2), mu=0, mu1=1, sigma=1, tau=2,
                   incomplete=IncMat2,verbose = T)
plot(c)
all.equal(a,c)

################################################################################
## plot.wlsPower #####
plot(wlsMixedPower(Cl=c(2,5),sigma=1,tau=0.1, mu0=0, mu1=1,
                            timepoints=5,design="parallel",verbose=T))


DesMat11 <- construct_DesMat(Cl=c(2,1))
CovMat11 <- construct_CovMat(3,3,1,2)
IncMat11 <- matrix(c(0,1,1,1,1,1,1,1,1),3,3)
IncMat22 <- matrix(c(1,1,1,1,0,1,1,1,1),3,3)
IncMat31 <- matrix(c(1,1,0,1,1,1,1,1,1),3,3)
IncMat2122 <- matrix(c(1,0,1,1,0,1,1,1,1),3,3)

wmp <- wlsMixedPower(DesMat=DesMat11, mu0=0, mu1=1, sigma=1, tau=2,
                          incomplete = NULL, verbose=TRUE)
wmpinc11 <- wlsMixedPower(DesMat=DesMat11, mu0=0, mu1=1, sigma=1, tau=2,
                          incomplete = IncMat11, verbose=TRUE)
wmpinc22 <- wlsMixedPower(DesMat=DesMat11, mu0=0, mu1=1, sigma=1, tau=2,
                          incomplete = IncMat22, verbose=TRUE)
wmpinc31 <- wlsMixedPower(DesMat=DesMat11, mu0=0, mu1=1, sigma=1, tau=2,
                          incomplete = IncMat31, verbose=TRUE)
wmpinc2122 <- wlsMixedPower(DesMat=DesMat11, mu0=0, mu1=1, sigma=1, tau=2,
                          incomplete = IncMat2122, verbose=TRUE)
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
## closed cohort

function(sigma, tau, eta, rho, alpha, K){

  tauNew <- sqrt(tau^2 + alph^2/K)

}


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


# ## wlsGlmmPower #####
sig <- .025
sig <- .03
tau <- 0.00262
wlsMixedPower(Cl=rep(2,6),mu0=0,mu1=.01,sigma=sqrt(sig*(1-sig)), tau=tau, N=250)
wlsMixedPower(Cl=rep(2,6),mu0=0,mu1=.01,sigma=sqrt(sig*(1-sig)/250), tau=tau, N=1)


debugonce(wlsGlmmPower)

mubar <- .025
sigSd <- sqrt(mubar*(1-mubar))
wlsGlmmPower (Cl=rep(2,6),mu0=.03,mu1=.02,            tau=.00262,N=250, marginal_mu = TRUE)
wlsGlmmPower (Cl=rep(2,6),mu0=.03,mu1=.02,            tau=.00262,N=250, marginal_mu = FALSE)
# wlsMixedPower(Cl=rep(2,6),EffSize = .01, sigma=sigSd, tau=.00262,N=250)
wlsMixedPower(Cl=rep(2,6),mu0=.03,mu1=.02,family="binomial", tau=.00262,N=250)

wlsGlmmPower (Cl=rep(2,6),mu0=.03,mu1=.02,            tau=.0262,N=250, marginal_mu = TRUE)
wlsGlmmPower (Cl=rep(2,6),mu0=.03,mu1=.02,            tau=.0262,N=250, marginal_mu = FALSE)
# wlsMixedPower(Cl=rep(2,6),EffSize = .01, sigma=sigSd, tau=.0262,N=250)

wlsGlmmPower (Cl=rep(2,6),mu0=.03,mu1=.02,            tau=.062,N=250, marginal_mu = TRUE)
wlsGlmmPower (Cl=rep(2,6),mu0=.03,mu1=.02,            tau=.062,N=250, marginal_mu = FALSE)
# wlsMixedPower(Cl=rep(2,6),EffSize = .01, sigma=sigSd, tau=.062,N=250)

wlsGlmmPower (Cl=rep(2,6),mu0=.03,mu1=.02,            tau=.262,N=250, marginal_mu = TRUE)
wlsGlmmPower (Cl=rep(2,6),mu0=.03,mu1=.02,            tau=.262,N=250, marginal_mu = FALSE)
# wlsMixedPower(Cl=rep(2,6),EffSize = .01, sigma=sigSd, tau=.262,N=250)

wlsGlmmPower (Cl=rep(2,6),mu0=.3,mu1=.2,                      tau=.62,N=20, marginal_mu = TRUE)
wlsGlmmPower (Cl=rep(2,6),mu0=.3,mu1=.2,                      tau=.62,N=20, marginal_mu = FALSE)
# wlsMixedPower(Cl=rep(2,6),EffSize = .1, sigma=sqrt(.25*.75),  tau=.62,N=20)


### compairing marginal-adjusted mu versus unadjusted (as in swCRTdesigns)
mu <- .235  ## mu conditional on tau=0
mu <- .285  ## marginal mu if muCond=.235
N  <- 60
sigSd <- sqrt(mu*(1-mu)/N)
tau   <- 3*sigSd
# wlsMixedPower(Cl=rep(2,6), EffSize=.05, sigma=sigSd, tau=tau, N=1)

N <- 1e6
binCond <- rbinom(N,1,p=inv_logit(logit(.2)+rnorm(N,0,.4)))
mean(binCond)

# wlsGlmmPower(Cl=c(1,1,1),mu0=0.04,mu1=0.02,tau=0.0,N=250)
# swPwr(swDsn(c(1,1,1)),mu0=.04,mu1=.02,tau=.0,eta=0,n=250,distn="binomial")
#
# wlsGlmmPower(Cl=rep(10,5),mu0=0.04,mu1=0.02,tau=0.01,N=1, verbose=F)
# swPwr(swDsn(rep(10,5)),mu0=.04,mu1=.02,tau=.01,eta=0,n=1,distn="binomial")
#
#
# ## binomial <-> gaussian analogy
# ## noch benoetigt fuer den "Wrapper"
# swPwr(swDsn(rep(25,4)),"binomial",200,0.03,0.025,tau=0.01,eta=0)
# sigtmp <- sqrt(.0275*.9725/200)
# microbenchmark::microbenchmark(
#   swPwr(swDsn(rep(25,4)),"gaussian",1,0.03,0.025,tau=0.01,eta=0,sigma=sigtmp)
#   ,
#   wlsMixedPower(EffSize=0.005,Cl=c(25,25,25,25),sigma=sigtmp,tau=0.01)[[1]]
#   ,times=100)
#
#

