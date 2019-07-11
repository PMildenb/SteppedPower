
## analytical power for iscb_poster

RUNS <- seq(0,1,length.out=11)
Power <- list()
for(tau in RUNS){
  Power[[paste0("_",tau)]] <- wlsGlmmPower(Cl=c(2,2,2,2,2,2),mu0=0.03, mu1=0.02, time_adjust="factor",
                                           tau_lin = tau,trt_delay=.5, N=250,verbose=T)$Power
}

Power


wlsGlmmPower(Cl=c(2,2,2,2,2,2),mu0=0.03, mu1=0.02, time_adjust="factor",
             tau_lin = 0.104,trt_delay=.5, N=250,verbose=T)$Power

