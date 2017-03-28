\name{subsetSDF}

\alias{subsetSDF}
\alias{subsetSDF.SpatialPointsDataFrame}
\alias{subsetSDF,SpatialPointsDataFrame-method}
\alias{subsetSDF.SpatialPixelsDataFrame}
\alias{subsetSDF,SpatialPixelsDataFrame-method}
\alias{subsetSDF.SpatialPolygonsDataFrame}
\alias{subsetSDF,SpatialPolygonsDataFrame-method}
\alias{subsetSDF.SpatialLinesDataFrame}
\alias{subsetSDF,SpatialLinesDataFrame-method}

\title{ Subsetting objects of class SpatialDataFrame}

\description{ Subsetting of objects of class \code{\link{SpatialDataFrame}}: selects according to the spatial objects or to the attribute(s), and returns an object of the same class. }

\usage{
subsetSDF(x, locations, data = names(x@data), ...)            
}
 
\arguments{
\item{x}{ object of class \code{SpatialDataFrame}}
\item{locations}{ integer vector: keep object at these spatial locations}
\item{data}{ integer vector or names of columns of the \code{x@data} to be kept }
\item{...}{further arguments to be passed to other methods}
}


\value{
An object of the same \code{SpatialDataFrame} class.
}

\note{
This is a wrapper function for the subsetting methods in all subclasses of \code{SpatialDataFrame}, with the parameters that are common to all classes; for the general subsetting methods refer to \code{\link[sp]{subset.Spatial}} for the classes imported from package \code{sp} and to \code{\link{subsetSDF.SpatialIndexDataFrame}} and \code{\link{subsetSDF.SpatialPolygridDataFrame}}.
}

\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
data(SIndexDF)
data(SPointsDF)
data(SPixelsDF)
data(SPolygridDF)
data(SPolygonsDF)
data(SLinesDF)

sub_Index = subsetSDF(SIndexDF, locations = c(1,3), data = "c")
sub_Points = subsetSDF(SPointsDF, locations = c(1,3), data = "z")
sub_Pixels = subsetSDF(SPixelsDF, locations = c(1,3), data = "z")
sub_Polygrid = subsetSDF(SPolygridDF, locations = c(1,3), data = "b")
sub_Polygons = subsetSDF(SPolygonsDF, locations = c(1,3), data = "a")
sub_Lines = subsetSDF(SLinesDF, locations = 1:2, data = "a")

}

%\keyword{}
