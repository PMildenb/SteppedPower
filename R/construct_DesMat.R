#'  construct_DesMat
#'
#' constructs the design matrix.
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD)
#' @param delay *not implemented*
#' @param timepoints numeric, scalar
#'
#' @return a matrix (for a stepped wedge design)
#' @examples construct_DesMat(Cl=c(2,0,1))
#'

construct_DesMat <- function(Cl,delay=NULL,design="SWD",timepoints=timepoints){

  SumCl         <- sum(Cl)
  if(design=="SWD"){
    sequences     <- length(Cl)
    timepoints    <- sequences + 1  ## modifikationen moeglich
    trt    <- matrix(0,sequences,timepoints) ; trt[upper.tri(trt)] <- 1
    trtBlk <- trt[rep(1:sequences,Cl),]
    trtvec <- as.numeric(t(trtBlk))
  }else if(design=="parallel"){
    if(length(Cl)!=2) {stop("Cannot handle length of vector Cl")}
    if(is.null(timepoints)){
      timepoints <- 1 ;  warning("timepoints unspecified. Defaults to 1.")}
    trtvec  <- c(rep(c(0,1),Cl*timepoints))
  }

  timeBlk <- suppressWarnings(cbind(1,rbind(0,diag(1,timepoints-1))))
  timeBlk <- timeBlk[rep(1:timepoints,SumCl),]
  DesMat  <- cbind(trtvec,timeBlk)
  return(DesMat)
}


