################################################################################
################################################################################
## speed tests
################################################################################
################################################################################

################################################################################
## Diagonal() vs diag() in construct_CovMat for gamma

library(Matrix)
gamma <- 11 ; timepoints <- 21 ; SumCl <- 40
M <- matrix(1, nrow=SumCl, ncol=SumCl)
microbenchmark::microbenchmark(
  (M %x% Diagonal(timepoints))*gamma,
  M %x% (Diagonal(timepoints)*gamma),
  M %x% Diagonal(timepoints)*gamma,
  kronecker(M,Diagonal(timepoints,gamma)),
  M %x% (diag(timepoints) *gamma),
  (M %x% diag(timepoints)) *gamma,
  times=5
)

################################################################################
##

dsnMat <- construct_DesMat(rep(10,50))

microbenchmark::microbenchmark(
  construct_DesMat(rep(10,50)),
  construct_DesMat(rep(10,50),time_adjust="linear"),
  construct_DesMat(rep(10,50),time_adjust="none")
)

################################################################################

dsnmatrix <- as.matrix(construct_DesMat(rep(1,70))$dsnmatrix)  ## dense matrix
DSNmatrix <- Matrix::Matrix(dsnmatrix,sparse=TRUE)             ## sparse matrix
CovMat <- construct_CovMat(SumCl=70,timepoints=71,sigma=1,tau=2)

microbenchmark::microbenchmark({
  tmpmat <- t(dsnmatrix) %*% Matrix::solve(CovMat)
  VarMat <- Matrix::solve(tmpmat %*% dsnmatrix)
},{
  tmpmat <- t(DSNmatrix) %*% Matrix::solve(CovMat)
  VarMat <- Matrix::solve(tmpmat %*% DSNmatrix)
},
  a <- compute_wlsPower(DesMat=construct_DesMat(rep(1,70)),
                   sigma=1, tau=2, EffSize=.1,verbose = TRUE)
,times=10)

################################################################################

microbenchmark::microbenchmark(
  wlsMixedPower(rep(1,3),EffSize=.1,sigma=1),
  wlsMixedPower(rep(1,3),EffSize=.1,sigma=1,tau=.1),
  wlsMixedPower(rep(1,3),EffSize=.1,sigma=1,tau=.1,eta=.2),
  wlsMixedPower(rep(1,3),EffSize=.1,sigma=1,tau=.1,eta=.2,rho=1)
)

microbenchmark::microbenchmark(
  wlsMixedPower(rep(2,10),EffSize=.1,sigma=1),
  wlsMixedPower(rep(2,10),EffSize=.1,sigma=1,tau=.1),
  wlsMixedPower(rep(2,10),EffSize=.1,sigma=1,tau=.1,eta=.2),
  wlsMixedPower(rep(2,10),EffSize=.1,sigma=1,tau=.1,eta=.2,rho=1)
)

microbenchmark::microbenchmark(
  wlsMixedPower(rep( 2,10),EffSize=.1,sigma=1,tau=.1,eta=.2,rho=1),
  wlsMixedPower(rep( 2,20),EffSize=.1,sigma=1,tau=.1,eta=.2,rho=1),
  wlsMixedPower(rep(10,10),EffSize=.1,sigma=1,tau=.1,eta=.2,rho=1),
  wlsMixedPower(rep(10,20),EffSize=.1,sigma=1,tau=.1,eta=.2,rho=1),
times=5)


################################################################################
## wlsMixedPower - CovMat/DesMat input #####
CM <- construct_CovMat(100,21,3,1)
DM <- construct_DesMat(rep(5,20))

microbenchmark::microbenchmark(
  construct_CovMat(100,21,3,1)
  ,
  wlsMixedPower(DesMat = DM, CovMat = CM, EffSize=5)
  ,
  wlsMixedPower(DesMat = DM, sigma=3, tau=1, EffSize=5)
  ,times=5)

## chol2inv instead of solve

microbenchmark::microbenchmark({
  t(DM$dsnmatrix)
},{
  a <- Matrix::solve(CM)
},{
  b <- Matrix::chol2inv(Matrix::chol(CM))
},{
#   c <- Rfast::spdinv(as.matrix(CM))
# },{
  tmpmat <- t(DM$dsnmatrix) %*% Matrix::solve(CM)
},{
  VarMat <- Matrix::solve(tmpmat %*% DM$dsnmatrix)
},unit="ms",times=20)

all.equal(a,b)

dim(a)
dim(b)
View(a)
View(b)
max(abs(b-a))



################################################################################
## really large designs

Cl_swd <- rep(5,60) ; EffSize <- .1 ; sigma <- 1 ; tau <- 1 ; design <- "SWD"
Cl_prl <- c(150,150) ; timepoints <- 61
system.time(
  pwrSWD <- wlsMixedPower(Cl=Cl_swd, design="SWD", EffSize=EffSize, sigma=sigma, tau=tau, verbose=TRUE))
system.time(
  pwrPRL <- wlsMixedPower(Cl=Cl_prl, design="parallel", EffSize=EffSize, sigma=sigma, tau=tau,
                          timepoints=timepoints, verbose=T))
dim(pwrPRL$DesignMatrix$dsnmatrix)
dim(pwrSWD$DesignMatrix$dsnmatrix)


system.time(
  HuHuPwrSWD <- swCRTdesign::swPwr(design=swCRTdesign::swDsn(rep(10,30)),distn="gaussian",
                                   n=1,mu0=0,mu1=.1,sigma=1,tau=1,eta=0,rho=0,gamma=0))

system.time(wlsMixedPower(Cl=rep(10,30), design="SWD", EffSize=EffSize, sigma=sigma, tau=tau, verbose=TRUE))



########################################################################################
## speed-test -> bei vielen Zeitp (mit wenigen Clustern) deutlich schneller.
Cl <- rep(1,45); sigtmp <- sqrt(.0275*.9725/200)

microbenchmark::microbenchmark(
  swCRTdesign::swPwr(design=swCRTdesign::swDsn(Cl),
                     distn="gaussian",
                     n=1, mu0=0.03, mu1=0.025,
                     tau=0.01, eta=0.01, rho=0, gamma=0, sigma=sigtmp)
  ,
  wlsMixedPower(EffSize=0.005,Cl=Cl,sigma=sigtmp,tau=0.01,eta=0.01)
  ,times=5)

## differenz  null
swCRTdesign::swPwr(design=swCRTdesign::swDsn(Cl),
                   distn="gaussian",
                   n=1, mu0=0.03, mu1=0.025,
                   tau=0.01, eta=0, rho=0, gamma=0, sigma=sigtmp) -  wlsMixedPower(EffSize=0.005,Cl=Cl,sigma=sigtmp,tau=0.01)[[1]]

## speed-test 2 viele Cluster, wenige Zeitpunkte
Cl <- rep(25,6);
microbenchmark::microbenchmark(
  swCRTdesign::swPwr(design=swCRTdesign::swDsn(Cl),
                     distn="gaussian",
                     n=1, mu0=0.03, mu1=0.025,
                     tau=0.01, eta=0, rho=0, gamma=0, sigma=sigtmp)
  ,
  wlsMixedPower(EffSize=0.005,Cl=Cl,sigma=sigtmp,tau=0.01)
  ,times=5)


########################################################################################
## contstruct toeplitz matrices efficently

library(Matrix)

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                    (1:n - 1))
  rho^exponent
}

n <- 800
microbenchmark::microbenchmark(
  ar1_cor(n,.8)
,
  stats::toeplitz(.8 ** c(0:(n-1)))
,
  Matrix::toeplitz(.8** c(0:(n-1)))
)

