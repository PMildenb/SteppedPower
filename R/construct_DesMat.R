#'  construct_DesMat
#'
#' constructs a
#'
#' @param Cl integer (vector), number of clusters per wave (in SWD)
#' @param delay *not implemented*
#'
#' @return a design matrix (for a stepped wedge design)
#' @examples construct_DesMat(Cl=c(2,0,1))
#'

construct_DesMat <- function(Cl,delay=NULL,design="SWD"){

  if(design=="SWD"){
    SumCl         <- sum(Cl)
    sequences     <- length(Cl)
    timepoints    <- sequences + 1  ## modifikationen moeglich
    trt    <- matrix(0,sequences,timepoints) ; trt[upper.tri(trt)] <- 1
    trtBlk <- trt[rep(1:sequences,Cl),]

    timeBlk <- cbind(1,rbind(0,diag(1,timepoints-1)))
    timeBlk <- timeBlk[rep(1:timepoints,SumCl),]

    DesMat  <- cbind(as.integer(t(trtBlk)),timeBlk)
  }else if(design=="SWD"){

  }

 return(DesMat)
}


