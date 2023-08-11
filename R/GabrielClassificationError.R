GabrielClassificationError = function(Data, ProjectedPoints, Cls, LC = NULL, 
                                      PlotIt = FALSE, Plotter = "native",
                                      Colors = NULL, LineColor= 'grey',
                                      main = "Name of Projection", mainSize = 24,
                                      xlab = "X", ylab = "Y", xlim, ylim,
                                      pch, lwd,
                                      Margin = list(t=50,r=0,l=0,b=0)){
  # GCE  = GabrielClassificationError(Data,ProjectedPoints,Cls) 
  # [GCE,GCEperPoint] = GabrielClassificationError(Data,ProjectedPoints,Cls) 
  # [GCE,GCEperPoint,nn,AnzNN,NNdists,HD] = GabrielClassificationError(OutputDelaunay,InputDistances,Cls) 
  # calculates GCE i.e the gabriel classification error
  #
  # INPUT
  # Data[1:n,1:d]					      Numeric matrix with n cases and d variables
  # ProjectedPoints[1:n,1:2]		Numeric matrix with 2D points in cartesian coordinates
  # Cls[1:n]						        Numeric vector with class labels
  # 
  # OPTIONAL
  # LC                        Numeric vector of two values determining grid size of the underlying projection
  # PlotIt                    Boolean: TRUE/FALSE => Plot/Do not plot (Default: FALSE)
  # Plotter                   Character with plot technique (native or plotly)
  # Colors                    Character vector of class colors for points
  # LineColor                 Character of line color used for edges of graph
  # main                      Character plot title
  # mainSize                  Numeric size of plot title
  # xlab                      Character name of x ax
  # ylab                      Character name of y ax
  # xlim                      Numeric vector with two values defining x ax range
  # ylim                      Numeric vector with two values defining y ax range
  # pch                       Numeric of point size (graphic parameter)
  # lwd                       Numeric of linewidth (graphic parameter)
  # Margin                    Margin of plotly plot
  #
  # OUPUT
  # GCE                            GabrielClassificationError 
  #                                NOTE the rest is just for development purposes
  # GCEperPoint(1:n)               unnormalized GCE of each point: GCE = mean( GCEperPoint)
  # nn                             the number of points in a relevant neghborhood: 0.5 * 85percentile(AnzNN)
  # AnzNN(1:n)                     the number of points with a gabriel graph neighborhood
  # NNdists(1:n,1:nn)              the distances within the relevant nehborhoot, 0 for inner cluster distances
  # HD(1:nn)                       HD = HarmonicDecay(nn) i.e weight function for the NNdists: GCEperPoint = HD*NNdists
  # IsInterDistance                Distances to the nn closest neighbors
  # gabrielDists                  Distance matrix implied by high dimensional
  #                                distances and the underlying gabriel
  #                                (Gabriel) graph
  # ProjectionGraphError           Plotly object in case, plotly is chosen
  # 
  # 
  # author: MT 07/2016, Sept 2016, Ausgabeparameteranpassung an ALU
  # Checked in 06/2018 by MT
  ## Catch Wrong Input
  if (!requireNamespace('DatabionicSwarm',quietly = TRUE)) {
    message(
      'Subordinate package (DatabionicSwarm) is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
    )
    return(
      "Subordinate package (DatabionicSwarm) is missing.
                Please install the package which is defined in 'Suggests'."
    )
  }
  
  if(!is.numeric(Data) | !is.matrix(Data)){
    stop("Data is not a numeric matrix")
  }
  
  if(!is.numeric(ProjectedPoints) | !is.matrix(ProjectedPoints)){
    stop("ProjectedPoints is not a numeric matrix")
  }
  
  if(!is.numeric(Cls) | !is.vector(Cls)){
    stop("Cls is not a numeric matrix")
  }
  
  if(is.list(ProjectedPoints)){
    stop('ProjectedPoints is a list not a matrix')
  }
  

  if(Plotter == "native"){
    if(missing(xlim)){
      xlim=c(min(ProjectedPoints[,1]),max(ProjectedPoints[,1]))
    }
    if(missing(ylim)){
      ylim=c(min(ProjectedPoints[,2]),max(ProjectedPoints[,2]))
    }
    if(missing(pch)){
      pch=20
    }
    if(missing(lwd)){
      lwd=0.1
    }
  }else if(Plotter == "plotly"){
    if(missing(xlim)){
      xlim=c(min(ProjectedPoints[,1])-0.5,max(ProjectedPoints[,1])+0.5)
    }
    if(missing(ylim)){
      ylim=c(min(ProjectedPoints[,2])-0.5,max(ProjectedPoints[,2])+0.5)
    }
    if(missing(pch)){
      pch=7
    }
    if(missing(lwd)){
      lwd=0.7
    }
  }
  
  c = ncol(ProjectedPoints)
  if (c > 3 | c < 2){
    stop(paste0('Wrong number of Columns of ProjectedPoints: ', c))
  }
  if(c == 3){ #Falls mit Key muss dieser in erster Spalte sein
    ProjectedPoints = ProjectedPoints[, 2:3]
  }
  AnzPunkte = length(Cls)
  if(length(Cls) < 3){
    stop('GabrielClassificationError: at least 3 points required')
  }
  
  #Calculate gabriel Graph in Outputspace  
  Gabriel=TRUE
  
  if(missing(LC) | is.null(LC)){ #Berechne gabriel Graphen, eigene Funktion, falls das mal auf CRAN soll
    del = Delaunay4Points(Points = ProjectedPoints, Gabriel = Gabriel, IsToroid = FALSE)
    if(!is.matrix(del)){
      del = del[[1]]
    }
  }else{
    # Im toroiden fall gibts die ESOM definition, wo x und y vertauscht sind
    #del = Delaunay4Points(Points = ProjectedPoints, IsToroid = TRUE, Grid = LC[c(2, 1)])
    del = Delaunay4Points(Points = ProjectedPoints, Gabriel = Gabriel, IsToroid = TRUE, LC = LC)
    if(!is.matrix(del)){
      del = del[[1]]
    }
  }
  
  #######################################################################
  HarmonicDecay = function(kj, n) {
    #  HD = HarmonicDecay(n)
    # calculates the Harmonic Decay Numbers HarmonicDecay(i,n) for i=1...n
    # HarmonicDecay(i,n) = sum(1/k) , k=i...n
    # Laws:
    # HarmonicDecay(1,n) = Harmonic number HN(n) = 1+1/2+...1/n
    #                      this approximates to:  log(n)+ 0.5772156649+ 1/(2*n)
    # HarmonicDecay(n,n) = 1/n
    #
    # INPUT
    # n      n > 0 natural number
    #
    # OUTPUT
    # HD(1:n)       HD(i) = HarmonicDecay(i,n)
    
    # MT Sept 2016
    
    if (n < 1)
      HD = NaN
    
    OneToN = c(1:n) # 1...n
    OneToN[OneToN > (kj-1)] = kj - 1
    
    HD = OneToN * NaN
    OneOverN = 1 / OneToN    # 1/1, ... 1/n
    for (i in 1:n) {
      HD[i] = sum(OneOverN[i:n])  # 1/i+ ...+ 1/n
    }
    #plot(OneToN,HD,ylab = 'HarmonicDecay(i,n)',xlab = 'i',main = 'HarmonicDecay')
    return(c(HD[1], HD[1:(length(HD)-1)]))
  }
  ####################################################################
  
  
  #Gewichte mit Euklid der gabriel kanten
  #Dist = as.matrix(dist(Data)) * del
  if(isSymmetric(Data)){
    Dist = Data * del
  }else{
    #storage.mode(del) = "numeric"
    Dist = as.matrix(dist(Data)) * del
  }
  ind             = which(Dist == 0, arr.ind = T)
  Dist[ind]       = NaN #Fuer spaeteres Sortieren sollten nicht direkt verbundene Punkte NaN sein, statt 0
  diag(Dist)      = NaN
  IsInterDistance = Dist * NaN     # IsInterDistance(i,j) =1 wenn die punkte i,j in verschiedene clustern
  gabrielDists   = Dist * NaN
  for(i in 1:nrow(Data)){
    SortInd              = order(Dist[, i], decreasing = F, na.last = T)
    currentClass         = Cls[i]
    cc                   = Cls[SortInd]
    IsInterDistance[, i] = currentClass != cc   # die inneren Distanzen Null setzen
    gabrielDists[, i]   = Dist[SortInd, i]
  }
  
  # jetzt die Anzahl der betrachten naesten Nachbarn und die Gewichtung ermitteln
  AnzNN = c()
  for(i in 1:nrow(Data)){
    AnzNN = c(AnzNN, sum(!is.nan(Dist[, i]))) #Anzahl echter Voronoi-Zelen
  }
  NN95 = round(stats::quantile(x = AnzNN, c(0.95)), 0)
  
  #nn      = round(NN95 / 2) # die Betrachtete nachbarschaft sind die Nachsten Nachbarn von 1...nn
  nn      = max(AnzNN)       # die Betrachtete nachbarschaft sind die Nachsten Nachbarn von 1...nn
  NNdists = IsInterDistance[1:nn, ]  # die ersten nn direkten Nachbahrn bzfl der Distanzengroesse
  
  for(i in 1:nrow(Data)){
    if(AnzNN[i] < nn){ #Setze 0, falls nicht direkt benachbahrt
      NotInterInd = seq(from = AnzNN[i] + 1, to = nn, by = 1)
      NNdists[NotInterInd, i] = 0
    }
  }
  
  # JM:
  # AnzNN, HarmonicDecay
  # N x max(nn)
  # lapply(AnzNN, HarmonicDecay)
  # 

  #tmpVar = lapply(AnzNN, HarmonicDecay)
  #tmpN = dim(ProjectedPoints)[1]
  #HarmonicMat = matrix(0, tmpN, nn)
  #
  #for(i in 1:length(tmpVar)){
  #  tmpVec1 = tmpVar[[i]]
  #  tmpVec2 = c(tmpVec1, rep(0, nn - length(tmpVec1)))
  #  HarmonicMat[i,] = tmpVec2
  #}
  
  # FIX: HarmonicDecay:

  NNdists     = t(NNdists)
  #HD         = HarmonicDecay(nn) # die Gewichtungsfuktion ausrechnen
  #GCEperPoint = NNdists %*% HD
  
  HD          = t(sapply(AnzNN, HarmonicDecay, nn))
  GCEpre      = NNdists * HD   # Gewicht * Distanzen: Vektorielle multiplikation
  GCEperPoint = apply(GCEpre, 1, sum)
  
  #GCEperPoint = apply(NNdists, 1, sum) #%*% HD   # Gewicht * Distanzen: Vektorielle multiplikation
  GCE         = mean(GCEperPoint, na.rm = T)  # GCE = GCEperPoint/AnzPunkte
  
  plotOut = "plotlyObject"
  if(isTRUE(PlotIt)){
  	 if(requireNamespace("DataVisualizations",quietly=TRUE)){
    if(Plotter == "native"){
      ind = which(GCEperPoint>0)
      DataVisualizations::PlotGraph2D(AdjacencyMatrix = del, Points = ProjectedPoints,
                                      lwd = lwd)
      points(ProjectedPoints[ind,],col="red")
    }else{
      #main=paste0(main, " - GCE:", round(GCE, 3))
      plotOut = DataVisualizations::PlotGraph2D(AdjacencyMatrix = del,
                                                Points          = ProjectedPoints,
                                                Cls             = Cls,
                                                Plotter         = "plotly",
                                                pch             = pch,
                                                lwd             = lwd,
                                                main            = main,
                                                mainSize        = mainSize,
                                                Colors          = Colors,
                                                LineColor       = LineColor,
                                                xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
      N = length(GCEperPoint)
      N=length(GCEperPoint)
      ErroxIdx  = which(GCEperPoint>0)
      NormalIdx = 1:N#setdiff(1:N, ErroxIdx)
      NoError   = round(GCEperPoint[NormalIdx], 2)
      Error     = round(GCEperPoint[ErroxIdx], 2)
	 if(requireNamespace("plotly",quietly=TRUE)){
      if(length(ErroxIdx) > 1){
        plotOut = plotly::add_markers(p = plotOut,
                                      x = ProjectedPoints[ErroxIdx,1],
                                      y = ProjectedPoints[ErroxIdx,2],
                                      text = Error, hoverinfo = 'text',
                                      type = "scatter",
                                      marker = list(size = pch, color = Colors[Cls[ErroxIdx]],
                                                    line = list(color = "red",
                                                                width = 2.5)))
      }
      plotOut = plotly::layout(p = plotOut, margin = Margin)
      plotOut = plotly::hide_legend(p = plotOut)
# "black",#
      print(plotOut)
	  }else{
	  message("GabrielClassificationError: for full plot please install plotly")
	  }
	  }
    }else{
	 message("GabrielClassificationError: for plot please install DataVisualizations")
	}
  }
  return(list("GCE"             = GCE,
              "GCEperPoint"     = GCEperPoint,
              "nn"              = nn,
              "AnzNN"           = AnzNN,
              "NNdists"         = NNdists,
              "HD"              = HD,
              "IsInterDistance" = IsInterDistance,
              "GabrielDists"   = gabrielDists,
              "ProjectionGraphError" = plotOut))
}    
           