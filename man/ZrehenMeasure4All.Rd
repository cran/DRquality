\name{ZrehenMeasure4All}
\alias{ZrehenMeasure4All}
\title{ A generalized version of the zrehen-measure which defines the neighbourhood by gabrielgraph and    is therefore not restricted to grid-based projections.}

\arguments{
\item{Data}{[1:n,1:d] points in input room with d attributes}
\item{Projection}{[1:n,1:2] projected points in output room, with index,x,y or index,line,column}
\item{width}{only necessary if toroid}
\item{height}{only necessary if toroid}
\item{isToroid}{are the points toroid?}
\item{isGrid}{is the grid a toroid?}

\item{plotGabriel}{plot the generated GabrielGraph}
}
\value{List with
\item{V$zrehen}{the raw zrehen measure}
\item{V$normedzrehen}{the zrehen measure normed by the number of neighbours}
\item{v$neighbourcounter}{the number of possible neighbours by which the zrehen measure is normed}
}
\description{
 A generalized version of the zrehen-measure which defines the neighbourhood by gabrielgraph and    is therefore not restricted to grid-based projections.}
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