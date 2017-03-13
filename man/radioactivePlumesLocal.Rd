\name{radioactivePlumes_local}

\alias{radioactivePlumes_local}
\docType{data}

\title{
Simulations of radioactive plumes around sources
}
\description{
The simulations show maps of pollution spreading from a single source in the centre. 
The data are the same as in the subsets A - C in \code{\link{radioactivePlumes_area}} but on a finer, telescopic resolution. The grid is fine (0.5km) close to the source and coarser (1km, 2km, 3km) farther away and reaches up to 100 km from the source.
}
\usage{
data(radioactivePlumes_local)
}
\format{
 \code{radioactivePlumes_local}: \code{\link{Simulations}}           
}
%\details{
%The \code{locations} contain \code{population} data (inhabitants per $km^2$) that was derived from the EEA population density maps (ref 1) by averaging over cells in the neibhbourhood of grid cells and \code{area} data ($km^2$). 
%}
%\references{
%  [1] European Environment Agency (2009).
%  \emph{Population density disaggregated with Corine land cover 2000}, 
%  \url{http://www.eea.europa.eu/data-and-maps/data/population-density-disaggregated-with-corine-land-cover-2000-2}
%}

\examples{
data(radioactivePlumes_local)

spplot(radioactivePlumes_local, zcol = 2)
}
\keyword{datasets}
