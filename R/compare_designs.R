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


compare_designs2 <- function(Designs, EffSize, family=gaussian(),
                             sigma, tau, eta=NULL, rho=NULL, sig.level=0.05){

}



# x <- wlsMixedPower(Cl=c(2,2,7,5,2), EffSize=.3, sigma=.5, tau=.5, eta=1.5, verbose=TRUE)
#
# plot_wls2 <- function(x){
#   if(!"WeightMatrix" %in% names(x)) stop("Please rerun wlsMixedPower with `verbose=TRUE` ")
#   wgt <- x$WeightMatrix
#   mx <- max(abs(wgt))
#   sumCl <- dim(wgt)[1]
#   timep <- dim(wgt)[2]
#   subp <- subplot(
#     plot_ly(data=data.frame(time=1:dim(wgt)[2], weight=colSums(abs(wgt))),
#             type="bar", x=~time, y=~weight, color=I("grey")) %>%
#       layout(xaxis=list(showticklabels=FALSE)),
#     plotly_empty(type="scatter",mode="marker"),
#     plot_ly(x=1:timep,y=1:sumCl,z=wgt,type="heatmap",
#             colors=colorRamp(c("steelblue","white","firebrick")),
#             xgap=.3,ygap=.3) %>%
#       colorbar(len=1,limits=c(-mx,mx)) %>% layout(yaxis = list(autorange = "reversed")),
#     plot_ly(data=data.frame(cluster=1:dim(wgt)[1], weight=rowSums(abs(wgt))),
#             type="bar", orientation="h",
#             y=~cluster, x=~weight, color=I("grey")) %>%
#       layout(yaxis=list(showticklabels=FALSE)),
#     nrows=2, heights=c(.2,.8), widths=c(.8,.2)
#   ) %>% layout(showlegend=FALSE)
#   subp
# }
# # debugonce(plot_wls2)
# plot_wls2(x)
#
# plot_wls3 <- function(x){
#   DesMat <- x$DesignMatrix
#   tauRange <- seq(.5*x$CovParams$tau, 2*x$CovParams$tau,length.out=5)
#   etaRange <- if(is.null(x$CovParams$eta)) 0:tauRange[2]
#
# }
#
