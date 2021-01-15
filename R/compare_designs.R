#' #' compare_designs
#' #'
#' #' a short wrapper that compares parallel, sw, and parallel+ baseline design
#' #'
#' #' @inheritParams wlsMixedPower
#' #'
#'
#' compare_designs <- function(Cl=NULL,
#'                             mu0,
#'                             mu1,
#'                             sigma,
#'                             tau,
#'                             family=gaussian(),
#'                             timepoints=NULL,
#'                             N=NULL,
#'                             sig.level=0.05,
#'                             DesMat=NULL,
#'                             trtDelay=NULL){
#'
#'   pwr_swd   <-  wlsMixedPower(mu0=mu0, mu1=mu1, sigma=sigma, tau=tau, timepoints=timepoints,
#'                             N=N, sig.level=sig.level,DesMat=DesMat, Cl=Cl, trtDelay=trtDelay,
#'                             design="SWD")
#'
#'   Cl_par <- c(ceiling(sum(Cl)/2),floor(sum(Cl)/2))
#'   if(is.null(timepoints)) timepoints <- length(Cl)+1
#'
#'   pwr_par   <-  wlsMixedPower(mu0=mu0, mu1=mu1, sigma=sigma, tau=tau, timepoints=timepoints,
#'                             N=N, sig.level=sig.level,DesMat=DesMat, Cl=Cl_par, trtDelay=trtDelay,
#'                             design="parallel")
#'   pwr_parBl <-  wlsMixedPower(mu0=mu0, mu1=mu1, sigma=sigma, tau=tau, timepoints=timepoints,
#'                             N=N, sig.level=sig.level,DesMat=DesMat, Cl=Cl_par, trtDelay=trtDelay,
#'                             design="parallel_baseline")
#'   return(data.frame("Stepped Wedge"= pwr_swd$Power,
#'                     "parallel"= pwr_par$Power,
#'                     "parallel+baseline"= pwr_parBl$Power,row.names="Power"))
#' }
#'
#' # Designs <- list(construct_DesMat(c(2,2,2)), construct_DesMat(c(3,0,3)))
#' # sigma <- 1 ; tau <- .5 ; eta <- .2 ; rho <- .1 ; EffSize <- .2
#' #' Title
#' #'
#' #' @inheritParams wlsMixedPower
#' #' @param Designs different designs to compare
#' #'
#' #' @return
#' #' @export
#'
#' compare_designs2 <- function(Designs, mu0,mu1, family=gaussian(),
#'                              sigma, tau, eta=NULL, rho=NULL, sig.level=0.05){
#'   out <- list()
#'   for(DesMat in Designs){
#'     out <- append(out, wlsMixedPower(DesMat=DesMat, mu0=mu0, mu1=mu1,
#'                          sigma=sigma, tau=tau, eta=eta, rho=rho, gamma))
#'   }
#'
#'   mapply(wlsMixedPower, DesMat=Designs,
#'                         mu0=mu0, mu1=mu1,
#'                         sigma=sigma, tau=tau, eta=eta, rho=rho, gamma)
#' }
#'
#'
#' ################################################################################
#'
