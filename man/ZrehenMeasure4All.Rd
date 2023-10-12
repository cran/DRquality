\name{ZrehenMeasure4All}
\alias{ZrehenMeasure4All}
\title{ZrehenMeasure4All}
\usage{ZrehenMeasure4All(Data, Projection, width, height, isToroid = FALSE,
isGrid = TRUE, plotGabriel = FALSE)}
\arguments{
\item{Data}{[1:n,1:d] points in input room with d attributes}
\item{Projection}{[1:n,1:2] projected points in output room, with index,x,y or index,line,column}
\item{width}{Numeric: only necessary if toroid}
\item{height}{Numeric: only necessary if toroid}
\item{isToroid}{Boolean: are the points toroid?}
\item{isGrid}{Boolean: is the grid a toroid?}
\item{plotGabriel}{Boolean: plot the generated GabrielGraph (TRUE) or not
(FALSE). Default: plotGabriel=FALSE.}
}
\value{List with
\item{V$zrehen}{the raw zrehen measure}
\item{V$normedzrehen}{the zrehen measure normed by the number of neighbours}
\item{v$neighbourcounter}{the number of possible neighbours by which the zrehen measure is normed}
}
\description{
A generalized version of the Zrehen-measure which defines the neighbourhood with
a Gabriel Graph and is therefore not restricted to grid-based projections.}
\author{Florian Lerch 07/2015}

\examples{
\donttest{
if(requireNamespace("FCPS")){
data(Hepta,package="FCPS")
projection=cmdscale(dist(Hepta$Data), k=2)
ZrehenMeasure4All(Hepta$Data,projection)$zrehen
}
}
\dontshow{
n=100
Data=cbind(rnorm(100),rnorm(100,2,3))
ZrehenMeasure4All(Data,cbind(runif(n),runif(n)))$zrehen
}
}