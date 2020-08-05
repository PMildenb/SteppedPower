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
#' @param trtmatrix
#' @param period
#'
#' @return a matrix (for a stepped wedge design)
#' @export
#'
#' @examples construct_DesMat(Cl=c(2,0,1))
#'

construct_DesMat <- function(Cl          =NULL,
                             trt_delay   =NULL,
                             design      ="SWD",
                             timepoints  =NULL,
                             time_adjust ="factor",
                             period      =NULL,
                             trtmatrix   =NULL,
                             timeBlk     =NULL){

  if(!is.null(trtmatrix)){
    trtMat <- trtmatrix
    if(inherits(trtMat,"matrix")){
      SumCl      <- nrow(trtMat)
      timepoints <- ncol(trtMat)
      Cl         <- table(do.call(paste,split(trtMat,col(trtMat))))  ## TODO: 1. add checks 2. find better alternative
    }else stop("trtmatrix must be a matrix. It is a ",class(trtMat))
  }else{
    trtMat  <- construct_trtMat(Cl            =Cl,
                                trt_delay     =trt_delay,
                                design        =design,
                                timepoints    =timepoints)
    timepoints <- dim(trtMat)[2]  ## trtMat has good heuristics for guessing timepoints (if not provided)
  }

  timeBlk <- construct_timeadjust(Cl          =Cl,
                                  timepoints  =timepoints,
                                  time_adjust =time_adjust,
                                  period      =period,
                                  timeBlk     =timeBlk)

  DesMat  <- list(dsnmatrix  = cbind(trt=as.numeric(t(trtMat)),timeBlk),
                  timepoints = timepoints,
                  SumCl      = sum(Cl),
                  trtMat     = trtMat)
  class(DesMat) <- append(class(DesMat),"DesMat")

  return(DesMat)
}

## Methods for class DesMat

#'  print.DesMat
#'
#' @param x d
#' @param ... Arguments to be passed to methods
#'
#' @method print DesMat
#' @export
#'
print.DesMat <- function(x, ...){
  cat("Timepoints         = ", x$timepoints,"\n")
  cat("Number of Clusters = ", x$SumCl,"\n")
  cat("Design matrix      = \n")
  print(x$dsnmatrix)
}


#' plot.DesMat
#'
#' @param x d
#' @param ... Arguments to be passed to methods
#'
#' @method plot DesMat
#' @export
#'

# x <- construct_DesMat(C=c(2,2,2,0,2,2,2),.5)
plot.DesMat <- function(x, ...){
  trt <- x$trtMat
  plot_ly(type="heatmap", colors=c("lightblue","red"),x=~1:dim(trt)[1],y=~1:dim(trt)[2],
          z=~trt, xgap=5, ygap=5) %>% layout(yaxis = list(autorange = "reversed"))
}



#' construct_trtMat
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
construct_trtMat <- function(Cl,trt_delay,design,timepoints=NULL){

  SumCl         <- sum(Cl)
  lenCl         <- length(Cl)

  if(design=="SWD"){
    if(is.null(timepoints)) timepoints <- length(Cl) + 1
    trt    <- matrix(0,lenCl,timepoints)
    trt[upper.tri(trt)] <- 1
    if(!is.null(trt_delay)){
      for(i in 1:length(trt_delay)){
        diag(trt[,-(1:i)]) <- trt_delay[i]  ## doesnt work if length(delay)>=length(timepoints)-1
      }
    }
  }else
  if(design=="parallel"){
    if(length(Cl)!=2) {stop("In construct_DesMat: Cl must be of length 2.")}
    if(is.null(timepoints)){
      if(is.null(trt_delay)){
        timepoints <- 1 ;  warning("timepoints unspecified. Defaults to 1.")
      }else{
        timepoints <- length(trt_delay)+1
        message("timepoints unspecified. Defaults to length(trt_delay)+1")
      }
    }
    trt     <- matrix(0,nrow=2,ncol=timepoints)
    trt[2,] <- c(trt_delay,rep(1,(timepoints-length(trt_delay))))
  }else
  if(design=="parallel_baseline"){
    if(length(Cl)!=2) {stop("In construct_DesMat: Cl must be of length 2.")}
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
    trt     <- matrix(0,nrow=2,ncol=timepoints)
    trt[2,] <- c(rep(0,timepoints01[1]),
                trt_delay,rep(1,(timepoints01[2]-length(trt_delay))))
  }
  trtMat <- as.matrix(trt[rep(1:lenCl,Cl),]) ## force matrix type (needed for timepoints==1)

  return(trtMat)
}


#' construct_timeadjust
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD)
#' @param timepoints numeric, scalar
#' @param time_adjust character, specifies adjustment for time periods. Defaults to "factor".
#' @param period number of timepoints per period. Defaults to `timepoints`
#'
#' @return What is returned? TODO
#' @export
#'
#' @examples
#'
construct_timeadjust <- function(Cl,timepoints=NULL,time_adjust=NULL,timeBlk=NULL,period=NULL){

  SumCl   <- sum(Cl)
  if(!is.null(timeBlk)) {
    timepoints <- dim(timeBlk)[1]
    timeBlks   <- timeBlk[rep(1:timepoints,SumCl),]
    return(timeBlks)
  }

  if(timepoints==1) time_adjust <- "none"
  if(time_adjust=="periodic" & is.null(period)) period <- timepoints

  timeBlk <- switch (time_adjust,
    factor   = cbind(1,rbind(0,Diagonal(timepoints-1)))[rep(1:timepoints,SumCl),],
    none     = matrix(rep(1,timepoints*SumCl)),
    linear   = matrix(rep(1:timepoints/timepoints,SumCl)),
    periodic = cbind(sin(0:(timepoints-1)*(2*pi/period)),
                     cos(0:(timepoints-1)*(2*pi/period)))[rep(1:timepoints,SumCl),]
  )

  return(timeBlks)
}
