#' compare_designs
#'
#' a short wrapper that compares parallel, sw, and parallel+ baseline design
#' @export

compare_designs <- function(EffSize,sigma,tau,family=gaussian(),timepoints=NULL,
                            N=NULL,sig.level=0.05,DesMat=NULL,Cl=NULL,trt_delay=NULL){

  pwr_swd   <-  wlsMixedPower(EffSize=EffSize, sigma=sigma, tau=tau, timepoints=timepoints,
                            N=N, sig.level=sig.level,DesMat=DesMat, Cl=Cl, trt_delay=trt_delay,
                            design="SWD")

  Cl_par <- c(ceiling(sum(Cl)/2),floor(sum(Cl)/2))
  if(is.null(timepoints)) timepoints <- length(Cl)+1

  pwr_par   <-  wlsMixedPower(EffSize=EffSize, sigma=sigma, tau=tau, timepoints=timepoints,
                            N=N, sig.level=sig.level,DesMat=DesMat, Cl=Cl_par, trt_delay=trt_delay,
                            design="parallel")
  pwr_parBl <-  wlsMixedPower(EffSize=EffSize, sigma=sigma, tau=tau, timepoints=timepoints,
                            N=N, sig.level=sig.level,DesMat=DesMat, Cl=Cl_par, trt_delay=trt_delay,
                            design="parallel_baseline")
  return(data.frame("Stepped Wedge"= pwr_swd$Power,
                    "parallel"= pwr_par$Power,
                    "parallel+baseline"= pwr_parBl$Power,row.names="Power"))
}
