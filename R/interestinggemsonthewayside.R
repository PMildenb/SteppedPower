## interesting observations


## groessere Cluster an den Rand bringt (kleine) Powerverbesserung (bekannt)
wlsMixedPower(EffSize=0.005,I=c(2,2,2,2),sigma=sigtmp,tau=0.01,N=rep(10,8))
wlsMixedPower(EffSize=0.005,I=c(2,2,2,2),sigma=sigtmp,tau=0.01,N=c(13,11,9,7,7,9,11,13))


## speed-test -> meine fkt ist deutlich schneller.
I <- rep(1,40); sigtmp <- sqrt(.0275*.9725/200)
microbenchmark(
  swPwr(swDsn(I),"gaussian",1,0.03,0.025,tau=0.01,eta=0,sigma=sigtmp)
  ,
  wlsMixedPower(EffSize=0.005,I=I,sigma=sigtmp,tau=0.01)
  ,times=10)
(swPwr(swDsn(I),"gaussian",1,0.03,0.025,tau=0.01,eta=0,sigma=sigtmp)
  -  wlsMixedPower(EffSize=0.005,I=I,sigma=sigtmp,tau=0.01))  ## differenz (nahezu) null

