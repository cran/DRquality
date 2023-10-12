\name{plotMeasureTandD}
\alias{plotMeasureTandD}
\title{Computes rank-based smoothed precision and recall}
\description{Compares the projection in pData with the original data in Data and
calculates trustworthiness and continuity of the projection for neighborhood
sizes ranging from 1 to the size of the neighborhood.}
\usage{plotMeasureTandD(TDmatrix, label = 'ProjectionMethod',
gPlotList = list(TW = ggplot2::ggplot(), DC = ggplot2::ggplot()), LineType = "solid",
Shape = 16, PointsPerE = 16)}
\arguments{
\item{TDmatrix}{Output of MeasureTundD() applied for a projection method.}
\item{label}{Title of plot.}
\item{gPlotList}{Settings for ggplot.}
\item{LineType}{Character - graphic parameter: Line type of ggplot.}
\item{Shape}{Integer: type of point}
\item{PointsPerE}{Numeric graphic parameter: Distance between markers on plot
line}
}
\value{ggplot object}
\author{
Michael Thrun
}
