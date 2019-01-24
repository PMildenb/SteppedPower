#'
#'  plot_wlsPower
#'
#' visual interpretation of cluster and period influences on intervention estimate
#' @param wlsPower Output of function wlsMixedPower
#'
#' @return three ggplot2 objects.
#' @export
#' @examples
#' # wlsPowerOut <- wlsMixedPower(EffSize=.1,I=c(1,2,1,0,1),sigma=.2,tau=0.05,N=10)
#' # plot_wlsPower(wlsPowerOut)
#'



plot_wlsPower <- function(wlsPower){
  WgtMat <- wlsPower$WeightMatrix
  HatData <- reshape2::melt(t(WgtMat[c(nrow(WgtMat):1),]))
  names(HatData) <- c("Time","Cluster","Weight")
  plotraw <- ggplot2::ggplot(HatData,ggplot2::aes_string("Time","Cluster")) +
    ggplot2::theme_minimal()  + ggplot2::scale_fill_gradient2(low="steelblue",mid="white",high="red") +
    ggplot2::geom_tile(ggplot2::aes_string(fill="Weight"),colour="white")


  Timeweights    <- colSums(abs(WgtMat))
  plotCluster <- ggplot2::ggplot(data.frame(Cluster=1:dim(WgtMat)[1],
                                            Clusterweights=rowSums(abs(WgtMat))),
                                 ggplot2::aes_string("Cluster","Clusterweights")) +
    ggplot2::geom_point() + ggplot2::theme_minimal()
  plotPeriods <- ggplot2::ggplot(data.frame(Periods=1:dim(WgtMat)[2],
                                             Weights=colSums(abs(WgtMat))),
                                 ggplot2::aes_string("Periods","Weights")) +
    ggplot2::geom_point() + ggplot2::theme_minimal()

  return(list(plotraw,plotCluster,plotPeriods))
}
