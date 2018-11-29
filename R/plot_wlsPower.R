#' 'plot_wlsPower'
#'
#'  visual interpretation of cluster and period influences on intervention estimate
#' @param wlsPower Output of function wlsMixedPower
#'
#' @return three ggplot2 objects.
#' @examples wlsPowerOut <- wlsMixedPower(EffSize=.1,I=c(1,2,1,0,1),sigma=.2,tau=0.05,N=10)
#' plot_wlsPower(wlsPowerOut)


plot_wlsPower <- function(wlsPower){
  HatMat <- wlsPower$HatMat
  HatData <- reshape2::melt(t(HatMat[c(nrow(HatMat):1),]))
  names(HatData) <- c("Time","Cluster","Weight")
  plotraw <- ggplot2::ggplot(HatData,aes(Time,Cluster)) +
    theme_minimal()  + scale_fill_gradient2(low="steelblue",mid="white",high="red") +
    geom_tile(aes(fill=Weight),colour="white")


  Timeweights    <- colSums(abs(HatMat))
  plotCluster <- ggplot2::ggplot(data.frame(Cluster=1:dim(HatMat)[1],
                                            Clusterweights=rowSums(abs(HatMat))),
                                 aes(Cluster,Clusterweights)) + geom_point() + theme_minimal()
  plotPeriods <- ggplot2::ggplot(data.frame(Periods=1:dim(HatMat)[2],
                                             Weights=colSums(abs(HatMat))),
                                 aes(Periods,Weights)) + geom_point() + theme_minimal()

  return(list(plotraw,plotCluster,plotPeriods))
}
