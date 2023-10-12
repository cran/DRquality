\name{MeasureTandD}
\alias{MeasureTandD}
\title{Trustworthiness and Discontinuity.}
\usage{MeasureTandD(Data, pData, NeighborhoodSize)}
\arguments{
\item{Data}{[1:n,1:d] points in input room with d attributes}
\item{pData}{[1:n,1:2] projected points in output room, with index,x,y or index,line,column}
\item{NeighborhoodSize}{Integer - sets the maximum number of neighbors to
calculate trustworthiness and continuity for.}
}
\value{
Numeric matrix [1:NeighborhoodSize,1:2] containing the trustworthiness values
in the first column and the discontinuity values in the second column.
}
\description{
In a trustworthy projection the visualized proximities hold in the original data
as well, whereas a continuous projection visualizes all proximities of the
original data.}
\author{Julian Märte}
\references{
Venna, J., & Kaski, S. (2005, September). Local multidimensional scaling with controlled tradeoff between trustworthiness and continuity. In Proceedings of 5th Workshop on Self-Organizing Maps (pp. 695-702).

Kaski, S., Nikkilä, J., Oja, M., Venna, J., Törönen, P., & Castrén, E. (2003). Trustworthiness and metrics in visualizing similarity of gene expression. BMC bioinformatics, 4(1), 1-13.
}
\examples{
\donttest{
if(requireNamespace("FCPS")){
data(Hepta,package="FCPS")
projection=cmdscale(dist(Hepta$Data), k=2)
MeasureTandD(Hepta$Data,projection, 2)
}
}
\dontshow{
n=100
Data=cbind(rnorm(100),rnorm(100,2,3))
MeasureTandD(Data,cbind(runif(n),runif(n)),2)
}
}

