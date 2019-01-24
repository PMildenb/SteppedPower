#'  construct_DesMat
#'
#' constructs the design matrix.
#'
#' Note: Unlike the usual notiation, the treatment is in the first column (for easier access by higher level functions).
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD)
#' @param delay *not implemented*
#' @param timepoints numeric, scalar
#'
#' @return a matrix (for a stepped wedge design)
#' @export
#'
#' @examples construct_DesMat(Cl=c(2,0,1))
#'

construct_DesMat <- function(Cl,delay=NULL,design="SWD",timepoints=timepoints){

  SumCl         <- sum(Cl)
  if(design=="SWD"){
    sequences     <- length(Cl)
    timepoints    <- sequences + 1  ## hier modifikationen moeglich
    trt    <- matrix(0,sequences,timepoints) ; trt[upper.tri(trt)] <- 1
    trtBlk <- trt[rep(1:sequences,Cl),]
    trtvec <- as.numeric(t(trtBlk))
  }else
  if(design=="parallel"){
    if(length(Cl)!=2) {stop("Cannot handle length of vector Cl")}
    if(is.null(timepoints)){
      timepoints <- 1 ;  warning("timepoints unspecified. Defaults to 1.")
    }
    trtvec  <- rep(c(0,1),Cl*timepoints)
  }else
  if(design=="parallel_baseline"){
    if(length(Cl)!=2) {stop("Cannot handle length of vector Cl")}
    if(length(timepoints)==1){
      timepoints01 <- c(1,timepoints-1)
      message(paste("assumes 1 baseline period and",timepoints-1,"parallel period(s)"))
    }else if(length(timepoints)==2){
      timepoints01 <- timepoints
      timepoints   <- sum(timepoints)
    }else if(is.null(timepoints)){
      timepoints01 <- c(1,1)
      timepoints   <- 2
      message("timepoints unspecified. Defaults to 1 baseline, 1 parallel period.")
    }
    trtvec  <- c(rep(0,timepoints*Cl[1]),
                 rep(rep(c(0,1),timepoints01),Cl[2]))
  }

  timeBlk <- suppressWarnings(cbind(1,rbind(0,diag(1,timepoints-1))))
  timeBlk <- timeBlk[rep(1:timepoints,SumCl),]
  DesMat  <- cbind(trtvec,timeBlk)
  return(DesMat)
}


