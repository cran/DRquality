\name{PrecisionAndRecall}
\alias{PrecisionAndRecall}
\title{Precision and Recall.}
\usage{PrecisionAndRecall(Data, pData, NeighborhoodSize = 20)}
\arguments{
\item{Data}{[1:n,1:d] points in input room with d attributes}
\item{pData}{[1:n,1:2] projected points in output room, with index,x,y or index,line,column}
\item{NeighborhoodSize}{Sets the 'effective number of neighbors' used to control
the width of the Gaussian, NeRV paper Seite 463 setzt Default auf 20}
}
\value{
Numeric matrix [1:NeighborhoodSize, 1:2] containing the precision values in the
first column and the recall values in the second column of the matrix.
}
\description{
Trade-off between missing similar points versus retrieving dissimilar points.}
\author{Felix Pape}
\examples{
\donttest{
if(requireNamespace("FCPS")){
data(Hepta,package="FCPS")
projection=cmdscale(dist(Hepta$Data), k=2)
PrecisionAndRecall(Hepta$Data,projection)
}
}
\dontshow{
n=100
Data=cbind(rnorm(100),rnorm(100,2,3))
PrecisionAndRecall(Data,cbind(runif(n),runif(n)))
}
}