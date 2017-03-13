\name{testDataTelescopic}

\alias{population_polygrid}
\alias{plumes_polygrid}

\docType{data}

\title{
Population and simulated plumes of gamma dose (rate) on a telescopic grid around a source of radiation.
}

\description{
Population data and simulated plumes of gamma dose.
The area covers 201 km times 201 km around a source. Values are given on a telescopic grid with resolution 0.5 km around the source and 1km, 2km, 3km farther away, it is given as \code{\link{SpatialPolygridDataFrame}}.

\code{population_polygrid}: population density

\code{plumes_polygrid}: total dose after plume passage


}
\usage{
data(population_polygrid)
data(plumes_polygrid)
}
\format{
 \code{population_polygrid}: \code{SpatialPolygridDataFrame}
 \code{plumes_polygrid}: \code{\link{SpatialPolygridDataFrame}}
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
  ~~ reference to a publication or URL from which the data were obtained ~~
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(plumes_polygrid)
spplot(plumes_polygrid)
}
\keyword{datasets}
