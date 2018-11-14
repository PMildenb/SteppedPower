#'  'construct_DesMat'
#'
#' @param I integer (vector), number of clusters per wave (in SWD)
#' @param delay not implemented
#'
#'
#'

construct_DesMat <- function(I,delay=NULL){
 SumCl         <- sum(I)
 sequences     <- length(I)
 timepoints    <- sequences + 1  ## modifikationen moeglich
 trt    <- matrix(0,sequences,timepoints) ; trt[upper.tri(trt)] <- 1
 trtBlk <- trt[rep(1:sequences,I),]

 timeBlk <- cbind(1,rbind(0,diag(1,timepoints-1)))
 timeBlk <- timeBlk[rep(1:timepoints,SumCl),]

 DesMat  <- cbind(as.integer(t(trtBlk)),timeBlk)

 return(DesMat)
}


