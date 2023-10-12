\name{SpearmansRho}
\alias{SpearmansRho}
\title{Calculates the error of a projection with spearman's rank correlation coefficient}
\usage{SpearmansRho(InputDists, OutputDists)}
\arguments{
\item{InputDists}{[1:d,1:d] numeric matrix with input distances}
\item{OutputDists}{[1:d,1:d] numeric matrix with output distances}
}
\value{rho}
\description{Calculates the error of a projection with spearman's rank
correlation coefficient}
\author{Julian Märte}
\examples{
\donttest{
if(requireNamespace("FCPS")){
data(Hepta,package="FCPS")
projection=cmdscale(dist(Hepta$Data), k=2)
SpearmansRho(as.matrix(dist(Hepta$Data)),as.matrix(dist(projection)))
}
}
\dontshow{
n=100
Data=cbind(rnorm(100),rnorm(100,2,3))
SpearmansRho(as.matrix(dist(Data)),as.matrix(dist(cbind(runif(n),runif(n)))))
}
}
