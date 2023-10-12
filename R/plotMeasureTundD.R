plotMeasureTandD = function(TDmatrix, label = 'ProjectionMethod',
                            gPlotList = list(TW = ggplot2::ggplot(), DC = ggplot2::ggplot()),
                            LineType="solid",Shape=16,PointsPerE=16) {
  #plotMeasureTundD(TDmatrix,label='ProjectionMethod',color='blue',gPlot=ggplot())
  #Plottet T und D von Venna/Kaski2001 als Kurvenverlauf ueber die ersten 50k
  # INPUT
  # TDmatrix[1:3,n]   		Output von MeasureTundD() eines Projektionsverfahrens
  # label								string, Bennenung der  Kurve
  # Optional
  # gPlotList         weiteren plotMeasureTundD uebergeben zum uewbereinander zeichnen
  # LineTypeq         linetyp als string
  # Shape             punktart als integer
  # PointsPerE        Abstand zwischen markierungspunkten auf linie, integer
  # Output
  #Ein ggPlot
  #author: MT 12/2015 und edited by MT 08/17
  #Example
  #   x=plotMeasureTundD(TundD[[1]],label=names[1])
  #   for(i in 2:length(names))
  #     x=plotMeasureTundD(TundD[[i]],label=names[i],gPlotList=x)

  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
    }else{ # Set up the page
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
      
      for (i in 1:numPlots) { # Make each plot, in the correct location
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                              layout.pos.col = matchidx$col))
      }
    }
  }
  
  knn = TDmatrix[, 1]#Siehe Readme in Doku
  

  ind=seq(from=1,to = length(knn),by=PointsPerE)
  
  df = data.frame(
    "KNN" = knn,
    "Trustworthiness" = (1 - TDmatrix[, 3]) * 100 ,
    "Projection" = factor(label)
  )#Siehe Readme in Doku vo NeRV
  #return(df)
  
  dfpoints= data.frame(
    "KNN" = knn[ind],
    "Trustworthiness" = (1 - TDmatrix[ind, 3]) * 100 ,
    "Projection" = factor(label)
  )#Siehe Readme in Doku vo NeRV
  #return(df)
  
  plt1 <-
    gPlotList$TW + ggplot2::geom_line(data = df,
                                      ggplot2::aes_string(x = "KNN",
                                                          y = "Trustworthiness",
                                                          colour = "Projection"),size=2,linetype=LineType) +
    ggplot2::geom_point(data = dfpoints,
                        ggplot2::aes_string(x = "KNN", y = "Trustworthiness", colour = "Projection"),
                        size = 5.5,
                        shape=Shape,
                        show.legend = FALSE,) +
    ggplot2::ylab('Trustworthiness in %')
  
  df2 = data.frame("KNN" = knn,
                   "Discontinuity" = (1 - TDmatrix[, 6]) * 100,
                   "Projection" = factor(label))#Siehe Readme in Doku TDmatrix
  dfpoints2 = data.frame("KNN" = knn[ind],
                         "Discontinuity" = (1 - TDmatrix[ind, 6]) * 100,
                         "Projection" = factor(label))#Siehe Readme in Doku TDmatrix
  
  plt2 = gPlotList$DC+ ggplot2::geom_line(data = df2,
                                          ggplot2::aes_string(x = "KNN",
                                                              y = "Discontinuity",
                                                              colour = "Projection"),size=2,linetype=LineType) +
    ggplot2::geom_point(data = dfpoints2,
                        ggplot2::aes_string(x = "KNN", y = "Discontinuity", colour = "Projection"),
                        size = 5.5, shape=Shape,
                        show.legend = FALSE,) +
    ggplot2::ylab('Discontinuity in %')
  
  
  multiplot(plt1, plt2, cols = 2)
  
  return(list(TW = plt1, DC = plt2))
}