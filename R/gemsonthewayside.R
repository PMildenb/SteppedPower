
########################################################################################
## interesting observations

########################################################################################
## groessere Cluster an den Rand bringt (kleine) Powerverbesserung (altbekannt)
sigtmp <- sqrt(.0275*.9725/200)
wlsMixedPower(EffSize=0.005,I=c(2,2,2,2),sigma=sigtmp,tau=0.01,N=rep(10,8))
wlsMixedPower(EffSize=0.005,I=c(2,2,2,2),sigma=sigtmp,tau=0.01,N=c(13,11,9,7,7,9,11,13))
wlsMixedPower(EffSize=0.005,I=c(2,2,2,2),sigma=sigtmp,tau=0.01,N=c(7,9,11,13,13,11,9,7))

########################################################################################
## speed-test -> meine fkt ist deutlich schneller.
library(microbenchmark) ; library(swCRTdesign)
I <- rep(1,40); sigtmp <- sqrt(.0275*.9725/200)

microbenchmark(
  swPwr(swDsn(I),"gaussian",1,0.03,0.025,tau=0.01,eta=0,sigma=sigtmp)
  ,
  wlsMixedPower(EffSize=0.005,I=I,sigma=sigtmp,tau=0.01)
  ,times=10)
(swPwr(swDsn(I),"gaussian",1,0.03,0.025,tau=0.01,eta=0,sigma=sigtmp)
  -  wlsMixedPower(EffSize=0.005,I=I,sigma=sigtmp,tau=0.01))  ## differenz (nahezu) null



########################################################################################
## Plot von gmds vortrag reloaded

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
wlsMixedPower(0.005,I=rep(25,4),sigma=sig.fct(200,0.025,0.025),tau=tau.fct(.4,.025))
wlsMixedPower(0.005,I=rep(25,4),sigma=sig.fct(200,0.03 ,0.03 ),tau=tau.fct(.4,.03 ))


## binomial wrapper fits nicely in between
## exp for tau=0
wlsMixedPower(0.005,I=rep(25,4),sigma=sig.fct(200,0.025,0.025),tau=0)
wlsMixedPower(0.005,I=rep(25,4),sigma=sig.fct(200,0.03 ,0.03 ),tau=0)
wlsGlmmPower(I=rep(25,4),mu0=0.03,mu1=0.025,tau=0,N=200)
swPwr(swDsn(rep(25,4)),"gaussian",1,0.03,0.025,tau=0,eta=0,sigma=sig.fct(200,0.025,0.025))
swPwr(swDsn(rep(25,4)),"gaussian",1,0.03,0.025,tau=0,eta=0,sigma=sig.fct(200,0.03,0.03))
swPwr(swDsn(rep(25,4)),"binomial",200,0.03,0.025,tau=0,eta=0) ## swPwr and wlsGlmm are ridiciously close

## for tau =!= 0
wlsMixedPower(0.005,I=rep(25,4),sigma=sig.fct(200,0.025,0.025),tau=tau.fct(.4,.025))
wlsMixedPower(0.005,I=rep(25,4),sigma=sig.fct(200,0.03 ,0.03 ),tau=tau.fct(.4,.03 ))
wlsGlmmPower(I=rep(25,4),mu0=0.03,mu1=0.025,tau=tau.fct(.4,0.025),N=200) ##
swPwr(swDsn(rep(25,4)),"binomial",200,0.03,0.025,tau=tau.fct(.4,.0275),eta=0)



########################################################################################
## What about package longpower ??
library(longpower)





















