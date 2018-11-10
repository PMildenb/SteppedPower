

## unorganized testing

## CorBlk
construct_CorBlk(timepoints=5, sigma=2, tau=2)

## CorMat
construct_CorMat(I=c(2,0,1), sigma=1, tau=0.3)

## DesMat
construct_DesMat(I=c(2,0,1))
construct_DesMat(I=c(2,2))

## wlsMixedPower
I <- c(2,2,2,2,2); sigma <- 4; tau <- 0.9 ; EffSize <- 1.5
wlsMixedPower(I=I,sigma=sigma,tau=tau,EffSize=EffSize)
swPwr(swDsn(I),distn="gaussian",1,0,EffSize,
      tau=tau,eta=0,rho=0,sigma=sigma)


wlsMixedPower(1,I=c(2,0,2,0,1),1,0.5)

