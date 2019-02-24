#'  construct_DesMat
#'
#' constructs the design matrix.
#'
#' Note: Unlike the usual notiation, the treatment is in the first column (for easier access by higher level functions).
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD)
#' @param trt_delay numeric (possibly vector), value(s) between 0 and 1 specifing the
#' intervention effect in the first (second ... ) intervention phase
#' @param design character, specifies the study design. Defaults to "SWD".
#' @param timepoints numeric, scalar
#' @param time_adjust character, specifies adjustment for time periods. Defaults to "factor".
#'
#' @return a matrix (for a stepped wedge design)
#' @export
#'
#' @examples construct_DesMat(Cl=c(2,0,1))
#'

construct_DesMat <- function(Cl, trt_delay=NULL, design="SWD", timepoints=NULL, time_adjust="factor"){

  trt_Lst    <- construct_trtvec(Cl=Cl, trt_delay=trt_delay, design=design, timepoints=timepoints)
  trtvec     <- trt_Lst[[1]]
  timepoints <- trt_Lst[[2]]

  timeBlk <- construct_timeadjust(Cl=Cl, timepoints=timepoints, time_adjust=time_adjust)

  DesMat  <- list(cbind(trtvec,timeBlk),timepoints,sum(Cl))
  names(DesMat) <- c("matrix","timepoints","SumCl")

  return(DesMat)
}


#' construct_trtvec
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD)
#' @param trt_delay numeric (possibly vector), value(s) between 0 and 1 specifing the
#' intervention effect in the first (second ... ) intervention phase
#' @param design character, specifies the study design. Defaults to "SWD".
#' @param timepoints numeric, scalar
#'
#' @return a vector trtvec
#' @export
#'
#' @examples
#'
#'
construct_trtvec <- function(Cl,trt_delay,design,timepoints){

  SumCl         <- sum(Cl)

  if(design=="SWD"){
    if(is.null(timepoints)) timepoints <- length(Cl) + 1
    trt    <- matrix(0,timepoints-1,timepoints)
    trt[upper.tri(trt)] <- 1
    if(!is.null(trt_delay)){
      for(i in 1:length(trt_delay)){
        diag(trt[,-(1:i)]) <- trt_delay[i]  ## doesnt work if length(delay)>=length(timepoints)-1
      }
    }
    trtBlk <- trt[rep(1:(timepoints-1),Cl),]
    trtvec <- as.numeric(t(trtBlk))
  }else
  if(design=="parallel"){
    if(length(Cl)!=2) {stop("In construct_DesMat: Cannot handle length of vector Cl")}
    if(is.null(timepoints)){
      if(is.null(trt_delay)){
        timepoints <- 1 ;  warning("timepoints unspecified. Defaults to 1.")
      }else{
        timepoints <- length(trt_delay)+1
        message("timepoints unspecified. Defaults to length(trt_delay)+1")
      }
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
      timepoints01 <- c(1,length(trt_delay)+1)
      timepoints   <- sum(timepoints01)
      message("timepoints unspecified. Defaults to 1 baseline,",length(trt_delay)+1, " parallel period(s).")
    }
    ctl_cluster <- rep(0,sum(timepoints01))
    trt_cluster <- c(rep(0,timepoints01[1]),
                     trt_delay,rep(1,(timepoints01[2]-length(trt_delay))))
    trtvec      <- c(rep(ctl_cluster,Cl[1]),rep(trt_cluster,Cl[2]))
  }


  # trtvec <- switch(design,
  #   SWD = {
  #     if(is.null(timepoints)) timepoints <- length(Cl) + 1
  #     trt    <- matrix(0,timepoints-1,timepoints)
  #     trt[upper.tri(trt)] <- 1
  #     if(!is.null(trt_delay)){
  #       for(i in 1:length(trt_delay)){
  #         diag(trt[,-(1:i)]) <- trt_delay[i]  ## doesnt work if length(delay)>=length(timepoints)-1
  #       }
  #     }
  #     trtBlk <- trt[rep(1:(timepoints-1),Cl),]
  #     as.numeric(t(trtBlk))
  #   },
  #   parallel = {
  #     if(length(Cl)!=2) {stop("In construct_DesMat: Cannot handle length of vector Cl")}
  #     if(is.null(timepoints)){
  #       if(is.null(trt_delay)){
  #         timepoints <- 1 ;  warning("timepoints unspecified. Defaults to 1.")
  #       }else{
  #         timepoints <- length(trt_delay)+1
  #         message("timepoints unspecified. Defaults to length(trt_delay)+1")
  #       }
  #     }
  #     ctl_cluster <- rep(0,timepoints)
  #     trt_cluster <- c(trt_delay,rep(1,(timepoints-length(trt_delay))))
  #     c(rep(ctl_cluster,Cl[1]),rep(trt_cluster,Cl[2]))
  #   },
  #   parallel_baseline = {
  #     if(length(Cl)!=2) {stop("In construct_DesMat: Cannot handle length of vector Cl")}
  #     if(length(timepoints)==1){
  #       timepoints01 <- c(1,timepoints-1)
  #       message(paste("assumes 1 baseline period and",timepoints-1,"parallel period(s)"))
  #     }else if(length(timepoints)==2){
  #       timepoints01 <- timepoints
  #       timepoints   <- sum(timepoints)
  #     }else if(is.null(timepoints)){
  #       timepoints01 <- c(1,length(trt_delay)+1)
  #       timepoints   <- sum(timepoints01)
  #       message("timepoints unspecified. Defaults to 1 baseline,",length(trt_delay)+1, " parallel period(s).")
  #     }
  #     ctl_cluster <- rep(0,sum(timepoints01))
  #     trt_cluster <- c(rep(0,timepoints01[1]),
  #                      trt_delay,rep(1,(timepoints01[2]-length(trt_delay))))
  #     c(rep(ctl_cluster,Cl[1]),rep(trt_cluster,Cl[2]))
  #   }
  # )

  return(list(trtvec,timepoints))
}


#' construct_timeadjust
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD)
#' @param timepoints numeric, scalar
#' @param time_adjust character, specifies adjustment for time periods. Defaults to "factor".
#'
#' @return
#' @export
#'
#' @examples
#'
construct_timeadjust <- function(Cl,timepoints,time_adjust){

  SumCl    <- sum(Cl)

  timeBlk <- switch (time_adjust,
    factor = cbind(1,rbind(0,diag(1,timepoints-1)))[rep(1:timepoints,SumCl),],
    none   = matrix(rep(1,timepoints*SumCl)),
    linear = matrix(rep(1:timepoints/timepoints,SumCl)),
    period = cbind(sin(0:(timepoints-1)*(2*pi/timepoints)),cos(0:(timepoints-1)*(2*pi/timepoints)))
  )

  return(timeBlk)
}
