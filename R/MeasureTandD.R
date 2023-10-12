MeasureTandD <- function(Data, pData, NeighborhoodSize){
  # MeasureTundD(Data, pData,NeighborhoodSize) 
  # Compares the projection in pData with the original data in Data
  # and calculates trustworthiness and continuity of the projection
  # for neighborhood sizes ranging from 1 to NeighborhoodSize
  #
  # INPUT
  # Data		          Matrix of original data
  # pData             Matrix of projected data
  # NeighborhoodSize  Sets the maximum number of neighbors to 
  #                   calculate trustworthiness and continuity for
  #
  # OUTPUT
  # Trustworthiness and Continuity for neighborhood sizes from 1 to NeighborhoodSize
  # (Actually each column contains the number 1 - x, where x is the measure. For example if the Output is 000593553, then
  #	the trustworthiness in this case is actually 1 - 0.000593553  = 0.9994.)
  #
  
  # Best-case and worst-case measures only come into play when you use a
  # projection method that tends to place many data points at exactly the same
  # spot (particularly the SOM). In that case a data point in the projection
  # will have several neighbors at exactly the same distance, and thus for some
  # k there will be more than one way to choose a set of k nearest neighbors.
  # `Best-case' trustworthiness is the value we get when we compute
  # trustworthiness with the nearest-neighbors set that gives the highest
  # trustworthiness; worst-case trustworthiness is the value we get when we use
  # the set that gives the lowest.

  if(!(is.matrix(Data) && is.matrix(pData) && length(NeighborhoodSize) == 1 && NeighborhoodSize %% 1 == 0) ){
    stop("Data and pData must be matrixes and NeighborhoodSize a natural number.")
  }
  if(dim(Data)[1] != dim(pData)[1]){
    stop("Data and pData have different numbers of rows.")
  }
  
  output = ProjectionBasedClustering::ContTrustMeasure(Data, pData, NeighborhoodSize)
  colnames(output) = c("Neighborhood size", "worst-case trustworthiness", "average trustworthiness", "best-case trustworthiness", "worst-case continuity", "average continuity", "best-case continuity")
  return(output)
}
