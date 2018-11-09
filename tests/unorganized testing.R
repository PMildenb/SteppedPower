

## unorganized testing


I <- c(2,2,2,2,2); sigma <- 4; tau <- 0.9 ; EffSize <- 1.5
d <- EffSize ; se <- sqrt(VarTrt)



swPwr(swDsn(c(2,2,2,2,2)),distn = "gaussian",1,0,1.5,0.9,0,0,sigma = 4)
View(swPwr)
