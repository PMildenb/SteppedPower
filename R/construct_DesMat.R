#'  construct_DesMat
#'
#' constructs the design matrix.
#'
#' Note: Unlike the usual notation, the treatment is in the first column (for easier access by higher level functions).
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD)
#' @param trtDelay numeric (possibly vector), value(s) between 0 and 1 specifying the
#' intervention effect in the first (second ... ) intervention phase
#' @param design character, specifies the study design. Defaults to "SWD".
#' @param timepoints numeric, scalar
#' @param timeAdjust character, specifies adjustment for time periods. Defaults to "factor".
#' @param period integer, only used if timeAdjust=="periodic". Defines the frequency of a periodic time trend.
#' @param trtmatrix an optional user defined matrix to define treatment allocation
#' @param timeBlk an optional user defined matrix that defines the time adjustment in one cluster.
#' Is repeated for every cluster.
#'
#' @return a matrix (for a stepped wedge design)
#' @export
#'
#' @examples construct_DesMat(Cl=c(2,0,1))
#'

construct_DesMat <- function(Cl          =NULL,
                             trtDelay   =NULL,
                             design      ="SWD",
                             timepoints  =NULL,
                             timeAdjust ="factor",
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
                                trtDelay     =trtDelay,
                                design        =design,
                                timepoints    =timepoints)
    timepoints <- dim(trtMat)[2]  ## trtMat has good heuristics for guessing timepoints (if not provided)
  }

  timeBlk <- construct_timeadjust(Cl          =Cl,
                                  timepoints  =timepoints,
                                  timeAdjust =timeAdjust,
                                  period      =period,
                                  timeBlk     =timeBlk)

  DesMat  <- list(dsnmatrix  = cbind(trt=as.numeric(t(trtMat)),timeBlk),
                  timepoints = timepoints,
                  Cl         = Cl,
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
  cat("Number of Clusters = ", sum(x$Cl),"\n")
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
  plot_ly(type="heatmap", colors=c("lightblue","red"),
          x=~1:dim(trt)[1],y=~1:dim(trt)[2],
          z=~trt, xgap=5, ygap=5) %>% layout(xaxis = list(title="Timepoints"),
                                             yaxis = list(title="Cluster",
                                                          autorange = "reversed"))
}



#' construct_trtMat
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD)
#' @param trtDelay numeric (possibly vector), value(s) between 0 and 1 specifying the
#' intervention effect in the first (second ... ) intervention phase
#' @param design character, specifies the study design. Defaults to "SWD".
#' @param timepoints numeric, scalar
#'
#' @return a matrix trtMat
#' @export
#'
#' @examples construct_trtMat(Cl=c(1,2,1), trtDelay=c(.2,.8), design="SWD")
#'
#'
construct_trtMat <- function(Cl,trtDelay,design,timepoints=NULL){

  SumCl         <- sum(Cl)
  lenCl         <- length(Cl)

  if(design=="SWD"){
    if(is.null(timepoints)) timepoints <- length(Cl) + 1
    trt    <- matrix(0,lenCl,timepoints)
    trt[upper.tri(trt)] <- 1
    if(!is.null(trtDelay)){
      for(i in 1:length(trtDelay)){
        diag(trt[,-(1:i)]) <- trtDelay[i]  ## doesnt work if length(delay)>=length(timepoints)-1
      }
    }
  }else if(design=="parallel"){
    if(length(Cl)!=2) {stop("In construct_DesMat: Cl must be of length 2.")}
    if(is.null(timepoints)){
      if(is.null(trtDelay)){
        timepoints <- 1 ;  warning("timepoints unspecified. Defaults to 1.")
      }else{
        timepoints <- length(trtDelay)+1
        message("timepoints unspecified. Defaults to length(trtDelay)+1")
      }
    }
    trt     <- matrix(0,nrow=2,ncol=timepoints)
    trt[2,] <- c(trtDelay,rep(1,(timepoints-length(trtDelay))))
  }else if(design=="parallel_baseline"){
    if(length(Cl)!=2) {stop("In construct_DesMat: Cl must be of length 2.")}
    if(length(timepoints)==1){
      timepoints01 <- c(1,timepoints-1)
      message(paste("assumes 1 baseline period and",timepoints-1,"parallel period(s).
                    If intended otherwise, argument timepoints must have length two."))
    }else if(length(timepoints)==2){
      timepoints01 <- timepoints
      timepoints   <- sum(timepoints)
    }else if(is.null(timepoints)){
      timepoints01 <- c(1,length(trtDelay)+1)
      timepoints   <- sum(timepoints01)
      message("timepoints unspecified. Defaults to 1 baseline,",length(trtDelay)+1, " parallel period(s).")
    }
    trt     <- matrix(0,nrow=2,ncol=timepoints)
    trt[2,] <- c(rep(0,timepoints01[1]),
                trtDelay,rep(1,(timepoints01[2]-length(trtDelay))))
  }else if (design=="crossover"){
    if(length(Cl)!=2) {stop("In construct_DesMat: Cl must be of length 2.")}
    if(length(timepoints)==1){
      timepoints01 <- c(floor(timepoints/2),ceiling(timepoints/2))
      message(paste("assumes", floor(timepoints/2) ,"AB period(s) and",
                    ceiling(timepoints/2),"BA period(s). If intended otherwise,
                    argument timepoints must have length two."))
    }else if (length(timepoints)==2){
      lenTp        <- timepoints
      timepoints   <- sum(timepoints)
    }else if(is.null(timepoints)){
      len          <- length(trtDelay)+1
      lenTp        <- c(len,len)
      timepoints   <- sum(timepoints01)
      message("timepoints unspecified. Defaults to", len, "AB period(s), and",
                len, " parallel period(s).")
    }
    trt   <- matrix(0, nrow=2, ncol=timepoints)
    vecAB <- c(trtDelay,rep(1,lenTp[1]-length(trtDelay)))
    vecBA <- c(trtDelay,rep(1,lenTp[2]-length(trtDelay)))

    trt[1,1:lenTp[1]]              <- vecAB
    trt[2,(lenTp[1]+1):timepoints] <- vecBA
  }
  trtMat <- as.matrix(trt[rep(1:lenCl,Cl),]) ## force matrix type (needed if timepoints==1)

  return(trtMat)
}


#' construct_timeadjust
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD)
#' @param timepoints numeric, scalar
#' @param timeAdjust character, specifies adjustment for time periods. Defaults to "factor".
#' @param period number of timepoints per period. Defaults to `timepoints` **experimental!**
#' @param timeBlk an optional user defined matrix that defines the time adjustment in one cluster.
#' Is repeated for every cluster.
#'
#' @return What is returned? TODO
#' @export

construct_timeadjust <- function(Cl,
                                 timepoints,
                                 timeAdjust="factor",
                                 period=NULL,
                                 timeBlk=NULL){

  SumCl   <- sum(Cl)
  if(!is.null(timeBlk)) {
    timepoints <- dim(timeBlk)[1]
    timeBlks   <- timeBlk[rep(1:timepoints,SumCl),]
    return(timeBlks)
  }

  if(timepoints==1) timeAdjust <- "none"
  if(timeAdjust=="periodic" & is.null(period)) period <- timepoints

  timeBlks <- switch (timeAdjust,
    factor   = cbind(1,rbind(0,diag(timepoints-1)))[rep(1:timepoints,SumCl),],
    none     = matrix(rep(1,timepoints*SumCl)),
    linear   = cbind(rep(1,timepoints*SumCl),rep(1:timepoints/timepoints,SumCl)),
    periodic = cbind(rep(1,timepoints),
                     sin(0:(timepoints-1)*(2*pi/period)),
                     cos(0:(timepoints-1)*(2*pi/period)))[rep(1:timepoints,SumCl),]
  )

  return(timeBlks)
}
