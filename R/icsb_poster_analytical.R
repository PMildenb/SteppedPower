
## analytical power for iscb_poster

RUNS <- seq(0,1,length.out=11)
Power <- list()
for(run in RUNS){
  Power[[paste0("_",run)]] <- wlsGlmmPower(Cl=c(2,2,2,2,2,2),mu0=0.03, mu1=.02, time_adjust="factor",
                                           tau_lin = run,trt_delay=.5, N=250,verbose=T)$Power
}
Power



RUNS <- seq(0.03,0.02,length.out=11)
Power <- list()
for(run in RUNS){
  Power[[paste0("_",run)]] <- wlsGlmmPower(Cl=c(2,2,2,2,2,2),mu0=0.03, mu1=run, time_adjust="factor",
                                           tau_lin = .1,trt_delay=.5, N=250,verbose=T)$Power
}
Power



