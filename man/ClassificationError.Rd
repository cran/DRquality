\name{ClassificationError}
\alias{ClassificationError}
\title{Classification Error (rate)}
\description{
Compares projected points to a given prior classification using knn classifier.
}
\usage{
ClassificationError(OutputDistances,Cls,k=5)
}
\arguments{
\item{OutputDistances}{[1:n,1:n] numeric matrix with distance matrix of
projected data.}
\item{Cls}{[1:n] Numeric vector containing class information.}
\item{k}{number of k nearest neighbors, in Venna 2010 set to 5 (here default)}
}
\details{
Projected points are evaluated by k-nearest neighbor classification accuracy (with k = 5), that is, each sample in the
visualization is classified by majority vote of its k nearest neighbors in the visualization, and the
classification is compared to the ground truth label. [Venna 2010].
}
\value{
List with three entries:
\item{Error}{Classification Error: 1-Accuracy[1]}
\item{Accuracy}{Accuracy }
\item{KNNCls}{[1:n]] cls of knn classifier}
}
\references{
Venna, J., Peltonen, J., Nybo, K., Aidos, H., and Kaski, S. Information retrieval perspective to nonlinear dimensionality reduction for data visualization. The Journal of Machine Learning Research, 11, 451-490.  (2010)

Gracia, A., Gonzalez, S., Robles, V., and Menasalvas, E. A methodology to compare Dimensionality Reduction algorithms in terms of loss of quality. Information Sciences, 270, 1-27.  (2014)
}
\author{
Michael Thrun
}
\note{
Here, the Outputdistances of the Projected points are used.
}
\examples{
\donttest{
if(requireNamespace("FCPS")){
data(Hepta,package="FCPS")
projection=cmdscale(dist(Hepta$Data), k=2)
ClassificationError(as.matrix(dist(projection)),Hepta$Cls)
}
}
\dontshow{
n=100
Data=cbind(rnorm(100),rnorm(100,2,3))
Cls=c(rep(1,50),rep(2,50))
projection=cmdscale(dist(Data), k=2)
ClassificationError(as.matrix(dist(projection)),Cls)
}
}
\keyword{Classification Error}% use one of  RShowDoc("KEYWORDS")
\keyword{ClassificationError}% __ONLY ONE__ keyword per line
\keyword{ClassificationErrorRate}% __ONLY ONE__ keyword per line
