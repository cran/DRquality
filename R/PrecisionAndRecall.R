PrecisionAndRecall <- function(Data, pData, NeighborhoodSize = 20){
  # PrecisionAndRecall(Data, pData, NeighborhoodSize)
  # Compares the projection in pData with the original data in Data
  # and calculates the smoothed recall and smooted precision.
  #
  # INPUT
  # Data		          Matrix of original data
  # pData             Matrix of projected data
  # NeighborhoodSize  Sets the 'effective number of neighbors' used
  # to control the width of the Gaussian, NeRV paper Seite 463 setzt Default auf 20
  #
  # OUTPUT
  # smoothed recall and smoothed precision
  #
  # AUTOR
  # FP
  
  if(!is.matrix(Data))
    stop("Data must be a matrix")
  if(!is.matrix(pData))
    stop("pData must be a matrix")
  if(NeighborhoodSize %% 1 != 0 || length(NeighborhoodSize) != 1 || NeighborhoodSize < 0)
    stop("NeighborhoodSize must be a single natural number")
  
  if(dim(Data)[1] == dim(pData)[1])
  
  return(ProjectionBasedClustering::KLMeasure(Data,pData,NeighborhoodSize))
  #return(c_klmeasure(Data,pData,NeighborhoodSize))
}
