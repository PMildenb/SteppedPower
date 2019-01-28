
########################################################################################
## interesting observations

########################################################################################
## groessere Cluster an den Rand bringt (kleine) Powerverbesserung (altbekannt)
sigtmp <- sqrt(.0275*.9725/200)
wlsMixedPower(EffSize=0.005,Cl=c(2,2,2,2),sigma=sigtmp,tau=0.01,N=rep(10,8))
wlsMixedPower(EffSize=0.005,Cl=c(2,2,2,2),sigma=sigtmp,tau=0.01,N=c(13,11,9,7,7,9,11,13))
wlsMixedPower(EffSize=0.005,Cl=c(2,2,2,2),sigma=sigtmp,tau=0.01,N=c(7,9,11,13,13,11,9,7))

wlsMixedPower(EffSize=0.001,Cl=rep(1,20),sigma=sigtmp,tau=01,
              N=c(1:10,10:1))
wlsMixedPower(EffSize=0.001,Cl=rep(1,20),sigma=sigtmp,tau=01,
              N=c(10:1,1:10))

########################################################################################
## speed-test -> bei vielen Zeitp (mit wenigen Clustern) deutlich schneller.
library(microbenchmark) ; library(swCRTdesign)
Cl <- rep(1,45); sigtmp <- sqrt(.0275*.9725/200)

microbenchmark(
  swPwr(swDsn(Cl),"gaussian",1,0.03,0.025,tau=0.01,eta=0,sigma=sigtmp)
  ,
  wlsMixedPower(EffSize=0.005,Cl=Cl,sigma=sigtmp,tau=0.01)
  ,times=5)
(swPwr(swDsn(Cl),"gaussian",1,0.03,0.025,tau=0.01,eta=0,sigma=sigtmp)
  -  wlsMixedPower(EffSize=0.005,Cl=Cl,sigma=sigtmp,tau=0.01)[[1]])  ## differenz (nahezu) null

########################################################################################
## difference in binom setting (for high diff mu0-mu1)

wlsGlmmPower(Cl=c(3,3,3,3),mu0=.5,mu1=.24,tau=0.1)
swPwr(swDsn(c(3,3,3,3)),mu0=.5,mu1=.24,tau=0.1,eta=0,n=1,distn="binomial")

########################################################################################
## KidSafe Setting Cl
swPwr(swDsn(c(2,2,2,2,2,2),extra.time=1,tx.effect=.5), distn="binomial",
      n=250, mu0=0.03, mu1=0.02, tau=0.00262, eta=0.00131, rho=-.25, retDATA=FALSE)


########################################################################################
## Plot von gmds vortrag reloaded
##  ich kann die simulationen noch nicht reproduzieren (!!)

tau.fct <- function(tau,mu){ tau*mu*(1-mu) }
sig.fct <- function(nInd,mu0,mu1){
  mu_bar <- mean(c(mu0,mu1))
  sigma  <- sqrt(mu_bar*(1-mu_bar)/nInd)
  return(sigma)}

swPwr(swDsn(rep(25,4)),"gaussian",1,0.03,0.025,tau=tau.fct(.515,.025),eta=0,sigma=sig.fct(200,0.025,0.025))
swPwr(swDsn(rep(25,4)),"gaussian",1,0.03,0.025,tau=tau.fct(.515,.03),eta=0,sigma=sig.fct(200,0.03,0.03))

# > rowMeans(p.25.4_0.515<0.05)[1:3]
# p.ztest.1a p.ztest.2b_min p.ztest.2b_max
# 0.8296         0.8236         0.7434
# > rowMeans(abs(p.25.4_0.515))[4:6]
# tx.est.1a tx.est.2b_min tx.est.2b_max
# 0.188413018   0.004973754   0.005005962


swPwr(swDsn(rep(25,4)),"gaussian",1,0.03,0.025,tau=tau.fct(.4,.025),eta=0,sigma=sig.fct(200,0.025,0.025))
swPwr(swDsn(rep(25,4)),"gaussian",1,0.03,0.025,tau=tau.fct(.4,.03),eta=0,sigma=sig.fct(200,0.03,0.03))

# > rowMeans(p.25.4_0.4<0.05)[1:3]
# p.ztest.1a p.ztest.2b_min p.ztest.2b_max
# 0.8258         0.8356         0.7704
# > rowMeans(abs(p.25.4_0.4))[4:6]
# tx.est.1a tx.est.2b_min tx.est.2b_max
# 0.186691687   0.005011320   0.004970041

## linear function can reproduce the linear results
wlsMixedPower(0.005,Cl=rep(25,4),sigma=sig.fct(200,0.025,0.025),tau=tau.fct(.4,.025))
wlsMixedPower(0.005,Cl=rep(25,4),sigma=sig.fct(200,0.03 ,0.03 ),tau=tau.fct(.4,.03 ))


## binomial wrapper fits nicely in between
## exp for tau=0
wlsMixedPower(0.005,Cl=rep(25,4),sigma=sig.fct(200,0.025,0.025),tau=0,verbose=F)
wlsMixedPower(0.005,Cl=rep(25,4),sigma=sig.fct(200,0.03 ,0.03 ),tau=0,verbose=F)
wlsGlmmPower(Cl=rep(25,4),mu0=0.03,mu1=0.025,tau=0,N=200,verbose=F)
swPwr(swDsn(rep(25,4)),"gaussian",1,0.03,0.025,tau=0,eta=0,sigma=sig.fct(200,0.025,0.025))
swPwr(swDsn(rep(25,4)),"gaussian",1,0.03,0.025,tau=0,eta=0,sigma=sig.fct(200,0.03,0.03))
swPwr(swDsn(rep(25,4)),"binomial",200,0.03,0.025,tau=0,eta=0) ## swPwr and wlsGlmm are ridiciously close

## for tau =!= 0
wlsMixedPower(0.005,Cl=rep(25,4),sigma=sig.fct(200,0.025,0.025),tau=tau.fct(.4,.025),verbose=F)
wlsMixedPower(0.005,Cl=rep(25,4),sigma=sig.fct(200,0.03 ,0.03 ),tau=tau.fct(.4,.03 ),verbose=F)
wlsGlmmPower(Cl=rep(25,4),mu0=0.03,mu1=0.025,tau=tau.fct(.4,0.025),N=200,verbose=F) ##
swPwr(swDsn(rep(25,4)),"binomial",200,0.03,0.025,tau=tau.fct(.4,.0275),eta=0)



########################################################################################
## What about package longpower ??
## -> different setting. not suitable


########################################################################################
## sampsize calculation works

devtools::install_github("PMildenb/SteppedPower")
library(SteppedPower) ; library(swCRTdesign)

## Wie viele Individuen pro Cluster braucht man bei 4 Clustern in 4 Sequenzen je nach Design fuer 90% pwr?
EffSize <- .5 ; sigma <- 1 ; tau <- .1
wlsMixedPower(Cl=rep(1,4), EffSize=EffSize, sigma=sigma, tau=tau, Power=.9, N_range=c(1,100))
wlsMixedPower(Cl=c(2,2), timepoints=5, design="parallel", EffSize=EffSize, sigma=sigma, tau=tau, Power=.9)

tau <- .15
wlsMixedPower(Cl=rep(1,4), EffSize=EffSize, sigma=sigma, tau=tau, Power=.9, N_range=c(1,100))
wlsMixedPower(Cl=c(2,2), timepoints=5, design="parallel", EffSize=EffSize, sigma=sigma, tau=tau, Power=.9)


## Wie viele Individuen pro Cluster braucht man bei 50 Clustern in 10 Sequenzen je nach Design?
wlsMixedPower(Cl=rep(5,10), EffSize=.1, sigma = 1, tau=.1, Power=.9, N_range=c(1,100))
wlsMixedPower(Cl=c(25,25), timepoints=11, design="parallel", EffSize=.1, sigma=1, tau=.1, Power=.9)
wlsMixedPower(Cl=c(25,25), timepoints=11, design="parallel_baseline", EffSize=.1, sigma=1, tau=.1, Power=.9)
wlsMixedPower(Cl=c(25,25), timepoints=c(2,3), design="parallel_baseline", EffSize=.1,
              sigma=1, tau=.1, Power=.9,N_range=c(1,100))



###########################################################################################
## compare different values for trt_delay

## parameters of setting to compare:
EffSize <- .4 ; sigma <- 1 ; tau <- .2 ; verbose <- FALSE

wlsMixedPower(Cl=rep(2,5),trt_delay=c(.3,.7),EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose,N_range=c(1,100))
wlsMixedPower(Cl=rep(2,5),trt_delay=c(.8),   EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose,N_range=c(1,100))
wlsMixedPower(Cl=rep(2,5),trt_delay=NULL,    EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose,N_range=c(1,100))

wlsMixedPower(Cl=c(5,5),timepoints=6,trt_delay=c(.3,.7),design="parallel",
              EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose)
wlsMixedPower(Cl=c(5,5),timepoints=6,trt_delay=c(.8),   design="parallel",
              EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose)
wlsMixedPower(Cl=c(5,5),timepoints=6,trt_delay=NULL,    design="parallel",
              EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose)



#################################
## Verteilungen hoeren ...

# switch()
  for(i in 1:1000){ beepr::beep() ; Sys.sleep(rexp(1,2)) }
  for(i in 1:1000){ beepr::beep() ; Sys.sleep(runif(1,0,1)) }

N <- 1e8
primes <- which(as.logical(matlab::isprime(1:N)))
prime_dist <- primes - c(0,primes[1:(length(primes)-1)])

car::densityPlot(prime_dist)

###########################################################################################
## special class for design matrix

setClass("design.matrix",contains="matrix",
         slots=c(SumCl="numeric"))
new("design.matrix",SumCl=4)

setClass("track", slots = c(x="numeric", y="numeric"))
myTrack <- new("track", x = -4:4, y = exp(-4:4))
slot(myTrack, "x")
slot(myTrack, "y") <- log(slot(myTrack, "y"))
utils::str(myTrack)











