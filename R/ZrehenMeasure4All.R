ZrehenMeasure4All <-
  function(Data,
           Projection,
           width,
           height,
           isToroid = FALSE,
           isGrid = TRUE,
           plotGabriel = FALSE) {
    #  Z=ZrehenMeasure4All(Data, Projection)
    # A generalized version of the zrehen-measure which defines the neighbourhood by gabrielgraph and
    #    is therefore not restricted to grid-based projections.
    # INPUT
    # Data[1:n,1:m]    n Projection in input room with m attributes
    # Projection[1:n,1:3]     n projected Projection in output room, with index,x,y or index,line,column
    # width               only necessary if toroid
    # height              only necessary if toroid
    # Optional
    # toroid              is the grid a toroid?
    # plotGabriel         plot the generated GabrielGraph
    # fast                use the faster but slightly incorrect implementation of gabrielgraph out of
    #			the spdep package
    # OUTPUT
    # list with
    # V$zrehen            the raw zrehen measure
    # V$normedzrehen      the zrehen measure normed by the number of neighbours
    # v$neighbourcounter  the number of possible neighbours by which the zrehen measure is normed
    # Author: FL 07/2015
    # 1.Editor: MT 01/2016
    # EXAMPLE
    # V = zrehenMeasure4All(Data,Projection)
    # zrehen = V$normedzrehen
    
    
    # Kritik/Anpassung am Zrehen Ma?
    # Die Normierung wird vorgeschlagen als die Anzahl aller Nachbarschaftspaare.
    # Da jedes Nachbarschaftspaar aber mehrere Intruder haben kann, sind Werte > 1 m?glich.
    # Au?erdem w?rden 2 Punkte die komplett au?einander projiziert werden durch unz?hlige
    # Intruder massiv ?berbewertet werden.
    # Hier z?hle ich daher pro Nachbarschaft maximal einen Intruder.
    
    fast = TRUE
    
    
    TileBM = function (BestMatches, Lines, Columns)
    {
      BestMatches = as.matrix(unname(BestMatches))
      if (ncol(as.matrix(BestMatches)) == 3) {
        BestMatches = cbind(
          as.numeric(BestMatches[, 1]),
          as.numeric(BestMatches[,
                                 2]),
          as.numeric(BestMatches[, 3])
        )
        TiledBestMatches <- rbind(
          BestMatches,
          cbind(BestMatches[,
                            1], (BestMatches[, 2] + Lines), BestMatches[, 3]),
          cbind(BestMatches[, 1:2], (BestMatches[, 3] + Columns)),
          cbind(
            BestMatches[, 1],
            (BestMatches[, 2] + Lines),
            (BestMatches[, 3] + Columns)
          )
        )
      }
      else {
        if (ncol(as.matrix(BestMatches)) == 2) {
          BestMatches = cbind(1:length(BestMatches[, 1]), BestMatches)
          BestMatches = cbind(
            as.numeric(BestMatches[, 1]),
            as.numeric(BestMatches[, 2]),
            as.numeric(BestMatches[,
                                   3])
          )
          TiledBestMatches <- rbind(
            BestMatches,
            cbind(BestMatches[,
                              1], (BestMatches[, 2] + Lines), BestMatches[,
                                                                          3]),
            cbind(BestMatches[, 1:2], (BestMatches[,
                                                   3] + Columns)),
            cbind(
              BestMatches[, 1],
              (BestMatches[,
                           2] + Lines),
              (BestMatches[, 3] + Columns)
            )
          )
        }
        else {
          stop("Error: Number of Rows is not 2 or 3, nothing could be done")
        }
      }
      return(TiledBestMatches = TiledBestMatches)
    }
    # define the euclidean norm
    norm_vec <- function(x)
      sqrt(sum(x ^ 2))
    
    buildAdjMatrix <- function(gabriel, size) {
      result <- matrix(0, nrow = size, ncol = size)
      
      x <- cbind(gabriel$from, gabriel$to)
      
      for (r in 1:nrow(x)) {
        pairs <- combn(x[r, ], 2)
        for (i in 1:ncol(pairs)) {
          result[pairs[1, i], pairs[2, i]] <- 1
          result[pairs[2, i], pairs[1, i]] <- 1
        }
      }
      
      result
    }
    
    neighboursByGabriel <-
      function(pointIndex, Projection, gabrielGraph) {
        # INPUT
        # pointIndex     Index of the point to which the neighbours should be returned
        # Projection[1:n,1:3]       n Projection with index,x,y
        # gabrielGraph   the calculated GabrielGraph
        # OUTPUT
        # the Projection out of "Projection" that are neighbours with pointIndex
        
        # vector isNeighbour equals 1 for every index with which the given point is neighboured
        isNeighbour <- gabrielGraph[pointIndex, ]
        
        # change 0 and 1 to FALSE and TRUE to use it as indices
        isNeighbour = as.logical(isNeighbour)
        
        matrix(Projection[isNeighbour, ], ncol = 3)
      }
    
    quadGridNeighbours <-
      function(pointIndex,
               Projection,
               width,
               height,
               toroid = TRUE) {
        # Example: quadGridNeighbours(c(40,30), 82, 50, TRUE)
        # NOTE: What should happen with Projection directly over each other?
        
        position = Projection[pointIndex, c(2, 3)]
        
        neighbours <- matrix(
          c(
            position - c(0, 1),
            position + c(0, 1),
            position - c(1, 0),
            position + c(1, 0),
            position
          ),
          byrow = TRUE,
          ncol = 2
        )
        
        # deal with Projection over border
        if (toroid == FALSE) {
          # only keep Projection within grid
          neighbours <- neighbours[which(neighbours[, 1] >= 1 &
                                           neighbours[, 1] <= height &
                                           neighbours[, 2] >= 1 &
                                           neighbours[, 2] <= width),]
        }
        else{
          # turn Projection to the other side through mod operator
          neighbours <- cbind(neighbours[, 1] %% height ,
                              neighbours[, 2] %% width)
          neighbours[which(neighbours[, 1] == 0), 1] = height
          neighbours[which(neighbours[, 2] == 0), 2] = width
        }
        
        # get the indices of Projection placed at neighbouring positions
        ind <- c(
          which(Projection [, 2] == neighbours[1, 1] &
                  Projection[, 3] == neighbours[1, 2]),
          which(Projection [, 2] == neighbours[2, 1] &
                  Projection[, 3] == neighbours[2, 2]),
          which(Projection [, 2] == neighbours[3, 1] &
                  Projection[, 3] == neighbours[3, 2]),
          which(Projection [, 2] == neighbours[4, 1] &
                  Projection[, 3] == neighbours[4, 2]),
          which(Projection [, 2] == neighbours[5, 1] &
                  Projection[, 3] == neighbours[5, 2])
        )
        
        # remove the point itself from its neighbourhood
        ind <- ind[which(ind != pointIndex)]
        matrix(Projection[ind, ], ncol = 3)
      }
    
    testQuadGridNeighbours <- function() {
      #Projection <- matrix(c(1,1,1,
      #                   2,1,2,
      #                   3,3,3,
      #                   4,1,1,
      #                   5,3,1), byrow=TRUE, ncol=3)
      Projection <- matrix(c(1, 1, 1), byrow = TRUE, ncol = 3)
      plot(Projection[, c(2, 3)])
      width = 3
      height = 3
      pointIndex = 1
      quadGridNeighbours(1, Projection, width, height, TRUE)
    }
    
    # ==== THE MEASURE STARTS HERE ====
    # Matrix that will contain the number of intruders for every pair of neighbours (if they are no neighbours then -1)
    intruder_m <- matrix(NaN, nrow(Projection), nrow(Projection))
    
    # add identifiers/indices to the Projection if not already given
    if (ncol(Projection) == 2)
      Projection <- cbind(1:nrow(Projection), Projection)
    
    # list to collect all intruders
    intruders = matrix(ncol = 3)
    
    # put Projection on a grid
    if (!isGrid) {
      if (isToroid)
        #dann mit bestmatches als voraussetzung
        tmpPoints <-
          DatabionicSwarm::ProjectedPoints2Grid(Projection, height, width)
      else
        #planar geht mit normalen punkten
        tmpPoints <- Projection
    } else
      tmpPoints <- Projection
    
    
    
    if (fast) {
      if (!requireNamespace('spdep', quietly = TRUE)) {
        message(
          'Subordinate package (spdep) is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
        )
        return(
          list(
            Object = "Subordinate package (spdep) is missing.
                Please install the package which is defined in 'Suggests'."
          )
        )
      }
      
      if (isToroid) {
        TiledBestMatches = TileBM(tmpPoints[], height, width)
        
        TiledX = TiledBestMatches[, 3]
        
        TiledY = height * 2 + 1 - TiledBestMatches[, 2]
        # So rechnet man BM koordinaten in XY um
        g <- spdep::gabrielneigh(TiledBestMatches[, c(2, 3)])
        
        TiledGabriel = buildAdjMatrix(g, nrow(Projection) * 4)
        
        Key = tmpPoints[, 1]
        
        AnzBestMatches = length(Key)
        
        Offset = c(0, 1, 2, 3) * AnzBestMatches
        gabrielGraph = TiledGabriel[c(1:AnzBestMatches), c(1:AnzBestMatches)]
        # initialisieren mit dem oberen viertel
        for (i in Offset) {
          for (j in Offset) {
            gabrielGraph = gabrielGraph + TiledGabriel[c(1:AnzBestMatches) + i, c(1:AnzBestMatches) +
                                                         j]
            
          }# end for j
        }# end for i
        # zurueckrechnen auf die ungekachelten werte
        #DelaunayAdjazenzMatrix <-UntileGraph4BM(TiledDelaunay, TiledBestMatches, Lines, Columns);
        gabrielGraph = (gabrielGraph > 0) * 1
        # alles auf 0/1 reduzieren
      }
      else{
        # no toroid
        g <- spdep::gabrielneigh(tmpPoints[, c(2, 3)])
        gabrielGraph = buildAdjMatrix(g, nrow(Projection))
      }
      if (plotGabriel) {
        AdjacencyMatrix = gabrielGraph
        Coordinates = Projection[, c(2, 3)]
        if (!requireNamespace('DataVisualizations', quietly = TRUE)) {
          message(
            'Subordinate package (DataVisualizations) is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
          )
          return(
            list(
              Object = "Subordinate package (DataVisualizations) is missing.
                Please install the package which is defined in 'Suggests'."
            )
          )
        }
        DataVisualizations::PlotGraph2D(AdjacencyMatrix, Coordinates, Plotter = "native")
      }
      } else{
        stop("ZrehenMeasure4All: Currently not available.")
        # not fast
        #print(tmpPoints)
        # if (isToroid)
        #   #dann mit bestmatches als voraussetzung
        #   gabrielGraph <-
        #     GraphAlgorithms::Gabriel4BestMatches(tmpPoints, c(height, width), isToroid, plotGabriel)$Gabriel
        # else
        #   #planar geht mit normalen punkten
        #   gabrielGraph = GraphAlgorithms::GabrielGraphMatrix(
        #     tmpPoints[, 2],
        #     tmpPoints[, 3],
        #     PlotIt = F,
        #     calcparallel = F,
        #     inst = 1
        #   )$Gabriel
      }
      
      # get the gabrielgraph to get neighbours from
      # if(toroid) gabrielGraph <- Gabriel4BestMatches(Projection, c(height,width), TRUE,plotGabriel)$Gabriel
      #if(!toroid){
      #  gabrielGraph <- GabrielGraphMatrix(Projection[,3],Projection[,2],plotGabriel)$Gabriel
      
      # remove the connection from a point to itself
      for (i in 1:nrow(gabrielGraph))
        gabrielGraph[i, i] = 0
      
      #}
      
      intruderCounter <- 0
      intruderCounter2 <- 1
      neighbourCounter <- 0
      intruderArray = array(NaN, dim = c(nrow(Projection), nrow(Projection), nrow(Projection)))
      
      for (i in 1:nrow(Projection)) {
        # iterate over every point
        # if toroid, get neighbours from the grid
        #if(toroid) neighbours <- quadGridNeighbours(i, Projection,width,height,TRUE)
        # get positions for every point that is neighboured through gabriel-graph
        #else neighbours <- neighboursByGabriel(i, Projection, gabrielGraph)
        
        neighbours <- neighboursByGabriel(i, Projection, gabrielGraph)
        
        if (nrow(neighbours) > 0)
          for (j in 1:nrow(neighbours)) {
            # get index of the neighbour (first element of the Projection returned by neighboursByGabriel
            neighbourIndex = neighbours[j, 1]
            if (i == neighbourIndex)
              break
            
            # break if this index is already covered so that every pair of neighbours is only counted once
            #if(neighbourIndex <= i) break;
            
            # count neighbour for normalization
            neighbourCounter = neighbourCounter + 1
            
            # Point in the middle between the straight line of point i and point neighbourIndex
            middlePoint = 0.5 * (Data[i, ] + Data[neighbourIndex, ])
            
            # distance of i to middlePoint / any other point must have a higher distance
            minDist = 0.5 * norm_vec(Data[i, ] - Data[neighbourIndex, ])
            
            # intruders for current pair of Projection
            tmpIntruder <- 0
            
            # find all "intruders"
            for (k in 1:nrow(Projection)) {
              # iterate over every point
              if (k == i ||
                  k == neighbours[j, 1])
                break
              # skip Projection i and j
              
              dist = norm_vec(Data[k, ] - middlePoint)
              
              # if intruder found, count it
              if (dist <= minDist) {
                tmpIntruder = tmpIntruder + 1
                #intruders = rbind(,intruders)
                intruderArray[i, neighbourIndex, k] = 1
              }
            }
            
            intruderCounter <- intruderCounter + tmpIntruder
            intruderCounter2 <-
              c(intruderCounter2, log(tmpIntruder + exp(1)))
            intruder_m[i, neighbourIndex] <- tmpIntruder
          }
      }
      
      # all_intruders / (all_observed_neighbourpairs * nr_of_possible_intruders_per_pair)
      # normedZrehen = intruderCounter/(neighbourCounter* (nrow(Projection)-2) ) # funktioniert nicht
      if (!missing(height))
        normedZrehen = intruderCounter / (height * width)
      else
        normedZrehen = intruderCounter / (min(Projection, na.rm = T) * max(Projection, na.rm =
                                                                             T))
      
      #MT Ergaenzung nach Diss
      ZrehenRight = intruder_m[lower.tri(intruder_m, diag = FALSE)]
      normedZrehen = sum(ZrehenRight, na.rm = T) / nrow(intruder_m) #Normierung durch Anzahl von Datenpunkten
      return(
        list(
          zrehen = intruderCounter,
          normedZrehen = normedZrehen,
          intruder_m = intruder_m,
          neighbourCounter = neighbourCounter,
          zre2 = intruderCounter2,
          intruderArray = intruderArray
        )
      )#, intruders = intruders) )
    }
    