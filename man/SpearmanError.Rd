\name{SpearmanError}
\alias{SpearmanError}
\title{Calculates the error of a projection with spearman's rank correlation coefficient.}
\arguments{
\item{VectorOfInputDists(1:n2)}{dissimilarities in Input Space between the n data points in vector form as produced by squareform(Dists(1:n,1:n))}
\item{VectorOfOutputDists(1:n2)}{dissimilarities in Output Space between the n data points in vector form as produced by squareform(Dists(1:n,1:n))}
}
\value{rho rank correlation coefficient}
\description{
Calculates the error of a projection with spearman's rank correlation
coefficient.}
\author{Florian Lerch}
\examples{
\donttest{
if(requireNamespace("FCPS")){
data(Hepta,package="FCPS")
projection=cmdscale(dist(Hepta$Data), k=2)
SpearmanError(as.matrix(dist(Hepta$Data)),as.matrix(dist(projection)))
}
}
\dontshow{
n=100
Data=cbind(rnorm(100),rnorm(100,2,3))
SpearmanError(as.matrix(dist(Data)),cbind(runif(n),runif(n)))
}
}

