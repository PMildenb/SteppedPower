
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
## difference in binom setting (for high diff mu0-mu1)

wlsGlmmPower(Cl=rep(1,3),mu0=.7,mu1=.3,tau=0.4,N=20,verbose=F)
swCRTdesign::swPwr(swCRTdesign::swDsn(cl=rep(1,3)),mu0=.7,mu1=.3,tau=0.4,eta=0,rho=0,n=20,distn="binomial")
## in extremen Situationen gibt es einige % Unterschied. swPwr-Methode ist konservativ.

########################################################################################
## KidSafe Setting Cl

devtools::install_github("PMildenb/SteppedPower")
library(SteppedPower)
wlsGlmmPower(Cl=c(2,2,2,0,2,2,2),mu0=0.03, mu1=0.02, trt_delay=.5, tau=0.00254,N=180,verbose=F)
wlsGlmmPower(Cl=c(2,2,2,1,1,2,2),mu0=0.03, mu1=0.02, trt_delay=.5, tau=0.00254,N=180,verbose=F)


swPwr(swDsn(c(2,2,2,2,2,2),extra.time=1), distn="binomial",
      n=250, mu0=0.03, mu1=0.02, tau=0.00262, eta=0.0, rho=0, retDATA=FALSE)
swPwr(swDsn(c(2,2,2,2,2,2),extra.time=1), distn="binomial",
      n=250, mu0=0.03, mu1=0.02, tau=0.00262, eta=0.0, rho=0, retDATA=FALSE)

tauKid <- 0.0044
wlsGlmmPower(Cl=c(6,6),timepoints=7,trt_delay=NULL,design="parallel", time_adjust="factor",
             mu0=0.03, mu1=0.02, tau=tauKid,N=250,verbose=F)
wlsGlmmPower(Cl=c(6,6),timepoints=1,trt_delay=NULL,design="parallel", time_adjust="factor",
             mu0=0.03, mu1=0.02, tau=tauKid,N=1750,verbose=F)

wlsGlmmPower(Cl=c(2,2,2,2,2,2),mu0=0.03, mu1=0.02, trt_delay=NULL, tau=0.00262,N=250,verbose=F)
wlsGlmmPower(Cl=c(3,3,0,0,3,3),mu0=0.03, mu1=0.02, trt_delay=NULL, tau=0.00262,N=250,verbose=F)
wlsGlmmPower(Cl=c(5,1,0,0,1,5),mu0=0.03, mu1=0.02, trt_delay=NULL, tau=0.00262,N=250,verbose=F)
wlsGlmmPower(Cl=c(6,0,0,0,0,6),mu0=0.03, mu1=0.02, trt_delay=NULL, tau=0.00262,N=250,verbose=T)

wlsGlmmPower(Cl=c(6,0,0,0,0,6),mu0=0.03, mu1=0.02, trt_delay=0, tau=0.00262,N=250,verbose=F)
wlsGlmmPower(Cl=c(6,6),timepoints=c(7),design="parallel_baseline",
             mu0=0.03, mu1=0.02, trt_delay=NULL, tau=0.00262,N=250,verbose=F)


sd <- sqrt( 0.00262^2 + 0.025*.975/250 )
pwr::pwr.t.test(d=0.01/sd,n=6)

sd <- sqrt( 0.00262^2 + 0.025*.975/1750 )
pwr::pwr.t.test(d=0.01/sd,n=6)


KidSafe     <- wlsGlmmPower(Cl=c(2,2,2,2,2,2),mu0=0.03, mu1=0.02, time_adjust="factor",
                        tau_lin = 1, trt_delay=.5, N=250,verbose=T)
KidSafe_lin <- wlsMixedPower(Cl=c(2,2,2,2,2,2),EffSize=.01,sigma=sqrt(.025*.975), time_adjust="factor",
                        tau=0.0254,trt_delay=.5, N=250,verbose=T)
KidSafe[[1]]
KidSafe_lin[[1]]

plot_wlsPower(KidSafe)[[1]]
plot_wlsPower(KidSafe_lin)[[1]]

########################################################################################
## KidSafe Setting 2019-03-06

## tatsaechliche (geschaetzte ) Aufnahmen aus Teiln.-praxen pro Quartal
N_cl <- 80

## Optionen zur FZP-optimierung
# - Inzidenz anheben
# - Perioden hinzufuegen
# - Zeiteffektsanpassung periodisch, o.ae.
# - staerkerer Effekt (moeglweise nur in Subgruppe)


## aktuelle FZP
wlsGlmmPower(Cl=c(2,2,2,2,2,2),mu0=0.03, mu1=0.02,
             tau=0.00262,trt_delay=.5, N=N_cl,verbose=F)

## hoehere Inzidenz (5% statt 3%)
wlsGlmmPower(Cl=c(2,2,2,2,2,2),mu0=0.05, mu1=2/3*0.05,
             tau=0.00262,trt_delay=.5, N=N_cl,verbose=F)

## hoehere Inzidenz (5% statt 3%)  +  eine Periode mehr
wlsGlmmPower(Cl=c(2,2,2,0,2,2,2),mu0=0.05, mu1=2/3*.05,
             tau=0.00262,trt_delay=.5, N=N_cl,verbose=F)

## hoehere Inzidenz (5% statt 3%)  +  zwei Perioden mehr
wlsGlmmPower(Cl=c(2,2,2,0,0,2,2,2),mu0=0.05, mu1=2/3*.05,
             tau=0.00262,trt_delay=.5, N=N_cl,verbose=F)

## hoehere Inzidenz (5% statt 3%)  +  eine Periode mehr  +
wlsGlmmPower(Cl=c(2,2,2,0,2,2,2),mu0=0.05, mu1=2/3*.05,
             tau=0.00262,trt_delay=.5, N=N_cl,verbose=F)


swCRTdesign::swPwr(swDsn(c(2,2,1,1,1,1,2,2),tx.effect=.5), distn="binomial",
      n=N_cl, mu0=0.03, mu1=0.02, tau=0.00262, eta=0.0, rho=0, retDATA=FALSE)

wlsGlmmPower(Cl=c(2,2,2,0,0,2,2,2),mu0=0.03, mu1=0.02,
                            tau=0.00262,trt_delay=.5, N=N_cl,verbose=F)
wlsGlmmPower(Cl=c(2,2,1,1,1,1,2,2),mu0=0.03, mu1=0.02,
             tau=0.00262,trt_delay=.5, N=N_cl,verbose=F)
wlsGlmmPower(Cl=c(2,2,2,0,2,0,2,2),mu0=0.03, mu1=0.02,
             tau=0.00262,trt_delay=.5, N=N_cl,verbose=F)


########################################################################################

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
EffSize <- .1 ; sigma <- 1 ; tau <- .03 ; verbose <- FALSE

wlsMixedPower(Cl=rep(2,5),trt_delay=c(.3,.7),EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose)
wlsMixedPower(Cl=rep(2,5),trt_delay=c(.8),   EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose)
wlsMixedPower(Cl=rep(2,5),trt_delay=NULL,    EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose)

wlsMixedPower(Cl=c(5,5),timepoints=6,trt_delay=c(.3,.7),design="parallel",
              EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose)
wlsMixedPower(Cl=c(5,5),timepoints=6,trt_delay=c(.8),   design="parallel",
              EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose)
wlsMixedPower(Cl=c(5,5),timepoints=6,trt_delay=NULL,    design="parallel",
              EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose)


###########################################################################################
## compare time as factor vs no time adjustment

## parameters of setting to compare:
EffSize <- .1 ; sigma <- 1 ; tau <- .03 ; verbose <- TRUE
ClSwd <- rep(2,8) ; ClPrl <- c(8,8) ; timepoints <- 9

wlsMixedPower(Cl=ClSwd,time_adjust="None",EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose)
wlsMixedPower(Cl=ClSwd,EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose)

wlsMixedPower(Cl=ClPrl,timepoints=timepoints,time_adjust="None",design="parallel",
              EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose)
wlsMixedPower(Cl=ClPrl,timepoints=timepoints,design="parallel",
              EffSize=EffSize,sigma=sigma,tau=tau,Power=.8,verbose=verbose)




###########################################################################################
## there are cases where one does not gain anything by adding clusters+timepoints in swd
## in terms of number of individual observations needed

N4 <- wlsMixedPower(Cl=rep(1,4),EffSize=.2,sigma=1,tau=.1,Power=.9)$N
N4*4*5

N5 <- wlsMixedPower(Cl=rep(1,5),EffSize=.1,sigma=1,tau=.1,Power=.9)$N
N5*5*6

N8 <- wlsMixedPower(Cl=rep(1,8),EffSize=.2,sigma=1,tau=.1,Power=.9)$N
N8*8*9

N4_3 <- wlsMixedPower(Cl=c(rep(1,4),0,0,0),EffSize=.2,sigma=1,tau=.1,Power=.9)$N
N4*4*8

wlsMixedPower(Cl=rep(1,8),EffSize=.1,sigma=1,tau=.1,N=N8,verbose=T)



###########################################################################################
## tau := .5 and not 1 to reach 90% power :)
## theory: quotient of taus is sqrt of quotient of timepoints

673/4
wlsMixedPower(Cl=c(673,673),timepoints=1,EffSize=.25,design="parallel",sigma=1,tau=1, N=1)
wlsMixedPower(Cl=c(169,169),timepoints=4,EffSize=.25,design="parallel",sigma=1,tau=.5,N=1)

SteppedPower::wlsMixedPower(Cl=c(6,15,3),EffSize=.25,sigma=(1/10),tau=1,N=1)
construct_DesMat(c(6,0))

###########################################################################################
##

SteppedPower::wlsMixedPower(Cl=c(2,2,2,0,2,2,2), EffSize=.01,
                            sigma=sqrt(.025*.975), tau=0.00254, eta=(0.0001),
                            trt_delay=.5, N=50)

swCRTdesign::swPwr(design=swCRTdesign::swDsn(c(2,2,2,0,2,2,2),.5), "gaussian",n=50,mu0=0,mu1=.01,
                   sigma=sqrt(.025*.975), tau=0.00254, eta=0.0001, rho=0, gamma=0)

###########################################################################################
##
TrtMatrix <- SteppedPower::construct_DesMat(rep(2,5))$trtMat
DESIGN    <- construct_DesMat(trtmatrix=TrtMatrix[rep(1:10,10),])
COVARIANCE<- construct_CovMat(10, 6, sigma=1, tau=.1, psi=5, INDIV_LVL=TRUE, N=10)
wlsPower(DesMat=DESIGN, CovMat = COVARIANCE, mu0=0, mu1=1)



