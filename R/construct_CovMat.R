


construct_CovBlk <- function(timepoints,sigma,tau){
  return(diag(sigma^2,timepoints) + tau^2)
}



construct_CovMat <- function(I,timepoints=NULL,sigma,tau,family=gaussian(),N=1){

 if(is.null(timepoints))  timepoints <- length(I)+1
 Sigmas <- rep(sigma,sum(I))

 CovBlks <- mapply(construct_CovBlk,sigma=Sigmas,
                   MoreArgs=list(timepoints=timepoints,tau=tau),SIMPLIFY = FALSE)
 return(Matrix::bdiag(CovBlks))
}


construct_CovBlk(5,1,0.2)
construct_CovMat(I=c(2,2),sigma=3,tau=0.7)


