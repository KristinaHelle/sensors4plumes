\name{testDataArtificial}

\alias{SIndexDF}
\alias{SPointsDF}
\alias{SPixelsDF}
\alias{SPolygridDF}
\alias{SPolygonsDF}
\alias{SLinesDF}

\docType{data}

\title{
Small artificial test data
}
\description{
Small artificial test data sets for demonstration, taken from the examples of the respective classes.
}
\usage{
data(SIndexDF)
data(SPointsDF)
data(SPixelsDF)
data(SPolygridDF)
data(SPolygonsDF)
data(SLinesDF)
}
\format{
 \code{SIndexDF}: \code{\link{SpatialIndexDataFrame}}           \cr
 \code{SPointsDF}: \code{\link[sp]{SpatialPointsDataFrame}}     \cr
 \code{SPixelsDF}: \code{\link[sp]{SpatialPixelsDataFrame}}     \cr
 \code{SPolygridDF}: \code{\link{SpatialPolygridDataFrame}}     \cr
 \code{SPolygonsDF}: \code{\link[sp]{SpatialPolygonsDataFrame}} \cr 
 \code{SLinesDF}: \code{\link[sp]{SpatialLinesDataFrame}}       \cr
}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
data(SPolygridDF)

spplot(SPolygridDF, zcol = "b", 
       sp.layout = list("sp.points", SpatialPoints(coordinates(SPolygridDF))))

data(SIndexDF)
spplot(SIndexDF)

data(SLinesDF)
areaSDF(SLinesDF)

}
\keyword{datasets}
