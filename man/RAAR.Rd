\name{RAAR}
\alias{RAAR}
\title{Rescaled average agreement rate}
\description{
Rescaled average agreement rate deduced by the co-ranking matrix from LCMC.
}
\usage{
RAAR(Data, ProjectedPoints, kmax = nrow(Data) - 2, PlotIt = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Data}{
Matrix containing n cases in rows, d variables in columns or a distance matrix which in this case has to be symmetric
}
\item{ProjectedPoints}{
 n by OutputDimension matrix containing coordinates of the Projection
}
\item{kmax}{
maximum of intervall 1:kmax of k nearest neighbors
}
\item{PlotIt}{
Optional: Should the output be plottet. Default: TRUE
}
}
\value{
A list containing:
\item{Raar}{Rescaled average agreement rate}
\item{Aar}{Average agreement rate}
}
\examples{
\donttest{
if(requireNamespace("FCPS")){
data(Hepta,package="FCPS")
projection=cmdscale(dist(Hepta$Data), k=2)
RAAR(Hepta$Data,projection,kmax=nrow(Hepta$Data)-2,PlotIt=TRUE)
}
}
\dontshow{
n=100
Data=cbind(rnorm(100),rnorm(100,2,3))
RAAR(Data,cbind(runif(n),runif(n)),kmax=nrow(Data)-2,PlotIt=TRUE)
}
}
\references{
Lee, J. A., Peluffo-Ordonez, D. H., & Verleysen, M. Multiscale stochastic neighbor embedding: Towards parameter-free dimensionality reduction. Paper presented at the Proceedings of 22st European Symposium on Artificial Neural Networks, Computational Intelligence And Machine Learning (ESANN) (2014).
}
\author{
Michael Thrun
}


