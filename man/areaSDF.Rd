\name{areaSDF}

\alias{areaSDF}
\alias{areaSDF.SpatialIndexDataFrame}
\alias{areaSDF.SpatialPointsDataFrame}
\alias{areaSDF.SpatialPixelsDataFrame}
\alias{areaSDF.SpatialPolygridDataFrame}
\alias{areaSDF.SpatialPolygonsDataFrame}
\alias{areaSDF.SpatialLinesDataFrame}
\alias{areaSDF,SpatialIndexDataFrame-method}
\alias{areaSDF,SpatialPointsDataFrame-method}
\alias{areaSDF,SpatialPixelsDataFrame-method}
\alias{areaSDF,SpatialPolygridDataFrame-method}
\alias{areaSDF,SpatialPolygonsDataFrame-method}
\alias{areaSDF,SpatialLinesDataFrame-method}

\title{
Areas of elements of SpatialDataFrame objects
}
\description{
Returns the areas of all elements in a SpatialDataFrame object, e.g. of all grid cells in a SpatialPixelsDataFrame or of all Polygons in a SpatialPolygonsDataFrame. If the elements are Points, 0s are returned with a warning.
}
\usage{
areaSDF(x)
}

\arguments{
  \item{x}{ \code{\link{SpatialDataFrame-class}}
}
}
\value{
Vector of the areas of elements of the \code{SpatialDataFrame} (same length as the \code{SpatialDataFrame}): \cr
\code{SpatialIndexDataFrame}, \code{SpatialPointsDataFrame}: 0s, with warning. \cr
\code{SpatialPixelsDataFrame}, \code{SpatialGridDataFrame}: cell size, repeated for each cell \cr
\code{SpatialPolygridDataFrame}: sum of cell sizes of all cells related to the respective index value \cr
\code{SpatialPolygonsDataFrame}: sum of area of all polygons related to the respective \code{Polygons} object \cr
\code{SpatialLinesDataFrame}: sum of length of lines related to the respective \code{Lines} object
}

\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
data(SIndexDF)
areaSDF(SIndexDF)

data(SPolygridDF)
areaSDF(SPolygridDF)

data(SLinesDF)
areaSDF(SLinesDF)
}

