\name{GabrielClassificationError}
\alias{GabrielClassificationError}
\title{
Gabriel Classification Error (GCE)
}
\description{
GCE searches for the k-nearest neighbors of the first gabriel neighbors
weighted by the Euclidean Distances of the Inputspace [Thrun et al, 2023]. GCE evaluates these
neighbors in the Output space. A low value indicates a better two-dimensional
projection of the high-dimensional Input space.
}
\usage{
GabrielClassificationError(Data,ProjectedPoints,Cls,LC,
PlotIt=FALSE,Plotter = "native", Colors = NULL,LineColor= 'grey',
main = "Name of Projection", mainSize = 24,xlab = "X", ylab = "Y", xlim, ylim,
pch,lwd,Margin=list(t=50,r=0,l=0,b=0))
}
\arguments{
\item{Data}{[1:n,1:d] Numeric matrix with n cases and d variables}
\item{ProjectedPoints}{[1:n,1:2] Numeric matrix with 2D points in cartesian coordinates}
\item{Cls}{[1:n] Numeric vector with class labels}
\item{LC}{Optional, Numeric vector of two values determining grid size of the underlying projection}
\item{PlotIt}{Optional, Boolean: TRUE/FALSE => Plot/Do not plot (Default: FALSE)}
\item{Plotter}{Optional, Character with plot technique (native or plotly)}
\item{Colors}{Optional, Character vector of class colors for points}
\item{LineColor}{Optional, Character of line color used for edges of graph}
\item{main}{Optional, Character plot title}
\item{mainSize}{Optional, Numeric size of plot title}
\item{xlab}{Optional, Character name of x ax}
\item{ylab}{Optional, Character name of y ax}
\item{xlim}{Optional, Numeric vector with two values defining x ax range}
\item{ylim}{Optional, Numeric vector with two values defining y ax range}
\item{pch}{Optional, Numeric of point size (graphic parameter)}
\item{lwd}{Optional, Numeric of linewidth (graphic parameter)}
\item{Margin}{Optional, Margin of plotly plot}
}
\details{
Gabriel Classification Error (GCE) makes an unbiased evaluation of distance- and
density-based structures which might be even non-linear separable. First, GCE
utilizes the information provided by a prior classification to assess projected
structures. Second, GCE applies the insights drawn from graph theory. Details
are described in [Thrun et al, 2023].}
\value{
list of several entries containing first the GCE itself as main result followed
by further entries which contain potential important information
\item{GCE}{Numeric: the 'Gabriel Classification Error'}
\item{GCEperPoint}{[1:n] unnormalized GCE of each point: GCE =
mean(GCEperPoint)}
\item{nn}{the number of points in a relevant neghborhood:
0.5 * 85percentile(AnzNN)}
\item{AnzNN}{[1:n] the number of points with a Gabriel graph neighborhood}
\item{NNdists}{[1:n,1:nn] the distances within the relevant neighborhood, 1 for
inter cluster distances and 0 for inner cluster distances}
\item{HD}{[1:nn] HD = HarmonicDecay(nn) i.e weight function for the NNdists:
GCEperPoint = HD*NNdists}
\item{IsInterDistance}{Distances to the nn closest neighbors.}
\item{GabrielDists}{Distance matrix implied by high dimensional distances and the underlying gabriel (Gabriel) graph }
\item{ProjectionGraphError}{Plotly object in case, plotly is chosen.}
}
\references{
[Thrun et al, 2023] Thrun, M.C, Märte, J., Stier, Q.: Analyzing Quality Measurements for Dimensionality Reduction, Machine Learning and Knowledge Extraction (MAKE), Vol 5., accepted, 2023.
}
\author{
Michael Thrun, Quirin Stier, Julian Märte
}

\examples{

if(requireNamespace("FCPS")){
data(Hepta,package="FCPS")
projection=cmdscale(dist(Hepta$Data), k=2)
GabrielClassificationError(Hepta$Data,projection,Hepta$Cls)$GCE
}

\donttest{
if(requireNamespace("FCPS")){
data(Hepta,package="FCPS")
projection=cmdscale(dist(Hepta$Data), k=2)
GabrielClassificationError(Hepta$Data,projection,Hepta$Cls)$GCE
}
}
\dontshow{
n=100
Data=cbind(rnorm(100),rnorm(100,2,3))
Cls=c(rep(1,50),rep(2,50))
GabrielClassificationError(Data,cbind(runif(n),runif(n)),Cls,LC=NULL)$GCE
}
}
\keyword{Classification Error}
\keyword{Gabriel Classification Error}
\keyword{GCE}
\keyword{Quality measure}
\keyword{QM}
\keyword{Quality measurement}
\keyword{Evaluation of projection methods}

