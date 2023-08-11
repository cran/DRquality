#Cmeasure <- function(Data,Projection,method = "pathlength",p=1,distance="euclidean"){
# # res <- Cmeasure(Data,Projection)
# # Calculate the C-Measure subtypes MinimalPathlength and MinimWiring
# # INPUT
# # Data    Vektor der Punkte in Eingaberaum
# # Projection    Vektor der Punkte in Ausgaberaum
# # k    Anzahl der Nachbarn in der Naehe
# # method   es wird nur Minimal Pathlength und Minimal Wiring implementiert
# # p     fuer Minimal Wiring 'wiring', Ausgaberaum, default=1
# # distance   Distanzmass in Nachbarnschaftberechnung
# # OUTPUT
# # MinimalPathlength    the calculated value
#    MinimWiring 			 the calculated value
# author MT

Cmeasure <- function(Data, Projection, k = 1) {

  if (!requireNamespace('igraph', quietly = TRUE)) {
    message(
      'Subordinate package (igraph) is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
    )
    return(
      list(
        Object = "Subordinate package (igraph) is missing.
                Please install the package which is defined in 'Suggests'."
      )
    )
  }
    if (!requireNamespace('cccd', quietly = TRUE)) {
    message(
      'Subordinate package (cccd) is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
    )
    return(
      list(
        Object = "Subordinate package (cccd) is missing.
                Please install the package which is defined in 'Suggests'."
      )
    )
  }
  
  KNNGraph= function (DistanceMatrix, k, Data) 
  {
    requireNamespace("cccd")
    requireNamespace("igraph")
    KNNGraphAdjMatrix = NULL
    tryCatch({
      if (missing(DistanceMatrix)) {
        result = cccd::nng(x = Data, k = k, mutual = TRUE, 
                           method = "Euclidean")
      }
      if (missing(Data)) {
        result = cccd::nng(dx = DistanceMatrix, k = k, mutual = TRUE, 
                           method = NULL)
      }
      KNNGraphAdjMatrix = igraph::get.adjacency(result, sparse = FALSE, 
                                                type = "both")
    }, error = function(e) {
      warning(paste0("KNNGraphAdjMatrix(): ", e))
      KNNGraphAdjMatrix = matrix(0, nrow(Data), ncol(Data))
    })
    return(KNNGraphAdjMatrix)
  }
  
  #k>1 nicht in papern definiert!
  #requireNamespace("Distances")
  #requireNamespace("GraphAlgorithms")
  InputD = as.matrix(dist(Data))
  OutputD = as.matrix(dist(Projection))
  spath = KNNGraph(OutputD, k = k)
  swiring = KNNGraph(InputD, k = k)
  return(c(
    MinimalPathlength = sum(InputD * spath),
    MinimWiring = sum(OutputD * swiring)
  ))
  
}
