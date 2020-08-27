#' compare_designs
#'
#' a short wrapper that compares parallel, sw, and parallel+ baseline design
#' @export

compare_designs <- function(EffSize,sigma,tau,family=gaussian(),timepoints=NULL,
                            N=NULL,sig.level=0.05,DesMat=NULL,Cl=NULL,trtDelay=NULL){

  pwr_swd   <-  wlsMixedPower(EffSize=EffSize, sigma=sigma, tau=tau, timepoints=timepoints,
                            N=N, sig.level=sig.level,DesMat=DesMat, Cl=Cl, trtDelay=trtDelay,
                            design="SWD")

  Cl_par <- c(ceiling(sum(Cl)/2),floor(sum(Cl)/2))
  if(is.null(timepoints)) timepoints <- length(Cl)+1

  pwr_par   <-  wlsMixedPower(EffSize=EffSize, sigma=sigma, tau=tau, timepoints=timepoints,
                            N=N, sig.level=sig.level,DesMat=DesMat, Cl=Cl_par, trtDelay=trtDelay,
                            design="parallel")
  pwr_parBl <-  wlsMixedPower(EffSize=EffSize, sigma=sigma, tau=tau, timepoints=timepoints,
                            N=N, sig.level=sig.level,DesMat=DesMat, Cl=Cl_par, trtDelay=trtDelay,
                            design="parallel_baseline")
  return(data.frame("Stepped Wedge"= pwr_swd$Power,
                    "parallel"= pwr_par$Power,
                    "parallel+baseline"= pwr_parBl$Power,row.names="Power"))
}

# Designs <- list(construct_DesMat(c(2,2,2)), construct_DesMat(c(3,0,3)))
# sigma <- 1 ; tau <- .5 ; eta <- .2 ; rho <- .1 ; EffSize <- .2
compare_designs2 <- function(Designs, EffSize, family=gaussian(),
                             sigma, tau, eta=NULL, rho=NULL, sig.level=0.05){
  out <- list()
  for(DesMat in Designs){
    out <- append(out, wlsMixedPower(DesMat=DesMat, EffSize=EffSize,
                         sigma=sigma, tau=tau, eta=eta, rho=rho, gamma))
  }

  mapply(wlsMixedPower, DesMat=Designs,
                        EffSize=EffSize,
                        sigma=sigma, tau=tau, eta=eta, rho=rho, gamma)
}



# x <- wlsMixedPower(Cl=c(2,2,2,2,1,1,2,2), trtDelay = .5, EffSize=.01, N=500,
#                    sigma=sqrt(0.0244), tau=0.005, eta=0.005, verbose=TRUE)
#
# # debugonce(plot_wls2)
# plot_wls2(x)
#
# plot_wls3 <- function(x){
#   p      <- x$CovParams
#   tauRange <- seq(0, 2*p$tau,length.out=5)
#   etaRange <- if(is.null(p$eta)) seq(0,max(tauRange),length.out=5) else seq(0,2*p$eta,length.out=5)
#
#   grid <- expand.grid(tau=tauRange,eta=etaRange)
#   pwr <- list()
#   for(i in 1:dim(grid)[1]){
#     pwr[[i]] <- wlsMixedPower(DesMat=x$DesignMatrix,EffSize=.01,sigma=p$sigma,N=200,
#                               tau=grid$tau[i], eta=grid$eta[i])[["Power"]]
#   }
#   tmp <- data.frame(grid,pwr=unlist(pwr))
#   plot_ly(type="scatter",mode="lines+marker",x=~tau, y=~pwr, split=~eta, data=tmp)
# }
#
#
