#'  construct_DesMat
#'
#' constructs the design matrix.
#'
#' Note: Unlike the usual notiation, the treatment is in the first column (for easier access by higher level functions).
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD)
#' @param trt_delay *implemented for swd, parallel and parallel w baseline*
#' @param timepoints numeric, scalar
#' @param time_adjust character, specifies adjustment for time periods. Defaults to "factor".
#'
#' @return a matrix (for a stepped wedge design)
#' @export
#'
#' @examples construct_DesMat(Cl=c(2,0,1))
#'

construct_DesMat <- function(Cl,trt_delay=NULL,design="SWD",timepoints=NULL,time_adjust="factor"){

  SumCl         <- sum(Cl)

  if(design=="SWD"){
    sequences     <- length(Cl)
    if(is.null(timepoints)) timepoints <- sequences + 1  ## hier modifikationen moeglich
    trt    <- matrix(0,sequences,timepoints) ; trt[upper.tri(trt)] <- 1
    if(!is.null(trt_delay)){
      for(i in 1:length(trt_delay)){
        diag(trt[,-(1:i)]) <- trt_delay[i]  ## doesnt work if length(delay)>=length(timepoints)-1
      }
    }
    trtBlk <- trt[rep(1:sequences,Cl),]
    trtvec <- as.numeric(t(trtBlk))
  }else
  if(design=="parallel"){
    if(length(Cl)!=2) {stop("In construct_DesMat: Cannot handle length of vector Cl")}
    if(is.null(timepoints)){
      if(is.null(trt_delay)){
        timepoints <- 1 ;  warning("timepoints unspecified. Defaults to 1.")
      }else
        timepoints <- length(trt_delay)+1 ;  message("timepoints unspecified. Defaults to length(trt_delay)+1")
    }
    ctl_cluster <- rep(0,timepoints)
    trt_cluster <- c(trt_delay,rep(1,(timepoints-length(trt_delay))))
    trtvec      <- c(rep(ctl_cluster,Cl[1]),rep(trt_cluster,Cl[2]))
  }else
  if(design=="parallel_baseline"){
    if(length(Cl)!=2) {stop("In construct_DesMat: Cannot handle length of vector Cl")}
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
    ctl_cluster <- rep(0,sum(timepoints01))
    trt_cluster <- c(rep(0,timepoints01[1]),
                     trt_delay,rep(1,(timepoints01[2]-length(trt_delay))))
    trtvec      <- c(rep(ctl_cluster,Cl[1]),rep(trt_cluster,Cl[2]))
  }

  if(time_adjust=="factor"){
    timeBlk <- suppressWarnings(cbind(1,rbind(0,diag(1,timepoints-1))))
    timeBlk <- timeBlk[rep(1:timepoints,SumCl),]
  } else
  if(time_adjust=="None"){
    timeBlk <- matrix(rep(1,timepoints*SumCl))
  }

  DesMat  <- list(cbind(trtvec,timeBlk),timepoints,SumCl)
  names(DesMat) <- c("matrix","timepoints","SumCl")
  return(DesMat)
}


