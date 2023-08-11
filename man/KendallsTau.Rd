\name{KendallsTau}
\alias{KendallsTau}
\title{
Statistical correlation by Kendall
}
\description{
Calculates the statistical correlation by Kendall. Basically a wrapper to
pcaPP::cor.fk.
}
\usage{
KendallsTau(InputDists, OutputDists)
}
\arguments{
  \item{InputDists}{
Matrix containing the distances of the first dataset.
}
  \item{OutputDists}{
Matrix containing the distances of the second dataset.
}
}
\value{
Equivalent to \code{\link[pcaPP]{cor.fk}}
}
\author{
Michael Thrun
}


\examples{

\donttest{
if(requireNamespace("FCPS")){
data(Hepta,package="FCPS")
InputDist=dist(Hepta$Data)
projection=cmdscale(InputDist, k=2)
KendallsTau(as.matrix(InputDist),as.matrix(dist(projection)))
}
}
\dontshow{
n=100
Data=cbind(rnorm(100),rnorm(100,2,3))
InputDist=dist(Data)
projection=cmdscale(InputDist, k=2)
KendallsTau(as.matrix(InputDist),as.matrix(dist(projection)))
}
}