\name{TopologicalCorrelation}
\alias{TopologicalCorrelation}
\title{Topological Correlation}
\usage{TopologicalCorrelation(Data,ProjectedPoints,type='norm',method,Kn=0)}
\arguments{
\item{Data}{[1:n, 1:d] a numeric matrix of the given n-dim. points: the rows represent the points and the columns represent the coordinates in the d-dim. space.}
\item{ProjectedPoints}{[1:n, 1:2] numeric matrix of Projected Points, if missing, method should be set!}
\item{method}{Determines whether the selected projections method for a given set of d-dim. points is a good choice. Therefor, a result of 1 means the selected projection method is good, and a result value of 0 means that the Visualization of the given Data in the two-dim. space doesnt fit the problem.}
\item{type}{How the paths in the adjacency matrix should be weighted. 'norm' representes path lenthgs of 1 and eucldidean represents the distance in the euclidean metric.}
\item{Kn}{k nearest neighbours in the graph. only needed in method is isomap and LocallyLinearEmbedding}
}
\value{TC value}
\description{
Calculates the Topological Correlation
}
\author{Hermann Tafo, Laukert Schlichting 07/2015}

\examples{
#requires DatabionicSwarm v2.2.1

\donttest{
if(requireNamespace("FCPS")){
#data(Hepta,package="FCPS")
#projection=cmdscale(dist(Hepta$Data), k=2)
#TopologicalCorrelation(Hepta$Data,projection)
}
}
\dontshow{
#n=100
#Data=cbind(rnorm(100),rnorm(100,2,3))
#TopologicalCorrelation(Data,cbind(runif(n),runif(n)))
}
}