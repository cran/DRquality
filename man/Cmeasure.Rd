\name{Cmeasure}
\alias{Cmeasure}
\title{C-Measure subtypes}
\description{Calculate the C-Measure subtypes of minimal path length and minimal
wiring}
\usage{Cmeasure(Data, Projection, k = 1)}
\arguments{
\item{Data}{[1:n,1:d] numerical matrix of points in input space.}
\item{Projection}{[1:n,1:2] numerical matrix of points in output space.}
\item{k}{Number of nearest neighbors, both measures set it always to k=1.}
}
\value{[1:2] Numerical vector of MinimalPathlength and MinimalWiring values.}
\author{Michael Thrun}

\examples{

\donttest{
if(requireNamespace("FCPS")){
data(Hepta,package="FCPS")
projection=cmdscale(dist(Hepta$Data), k=2)
Cmeasure(Hepta$Data,projection)
}
}
\dontshow{
n=100
Data=cbind(rnorm(100),rnorm(100,2,3))
Cmeasure(Data,cbind(runif(n),runif(n)))
}
}