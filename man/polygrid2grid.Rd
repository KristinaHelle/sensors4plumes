\name{polygrid2grid}

\alias{polygrid2grid}

\title{Coerce SpatialPolygridDataFrame into SpatialGridDataFrame and geoTiff file}

\description{
Coerces \code{\link{SpatialPolygridDataFrame}} into \code{\link{SpatialGridDataFrame}} by copying the data to all associated grid cells. It may write the result to a geoTiff file.
}
\usage{
polygrid2grid(obj, zcol = NA, returnSGDF = TRUE, geoTiffPath)
}
\arguments{
  \item{obj}{\code{SpatialPolygridDataFrame-class}
}
  \item{zcol}{names or numbers of the data columns to be used
}
  \item{returnSGDF}{if the generated SpatialGridDataFrame is to be returned
}
  \item{geoTiffPath}{filename where to save the result as geoTiff (without suffix, '.tif' is added automatically); if missing, no geoTiff is generated
}
}
\details{
The main reason to turn a \code{SpatialPolygridDataFrame} into a \code{SpatialGridDataFrame} is for plotting. It can make sense to generate a (multilayer) geoTiff instead of returning the transformed data into the workspace.
}
\value{
A \code{SpatialGridDataFrame} with the same values assigned as in the \code{SpatialPolygridDataFrame} by copying the data to all associated grid cells. If \code{returnSGDF = FALSE} it returns \code{TRUE}.
}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
data(SPolygridDF)

# return SpatialGridDataFrame
SGridDF1 = polygrid2grid(SPolygridDF, zcol = "b")

# generate geoTiff 
polygrid2grid(SPolygridDF, returnSGDF = FALSE, 
              geoTiffPath = "SPolygridDF1")

# plot SpatialGridDataFrame
spplot(SGridDF1)

# view geoTiff with functions from 'raster' and delete it
SGridDF2 = brick("SPolygridDF1.tif")
plot(SGridDF2)
rm(SGridDF2)
file.remove("SPolygridDF1.tif")
}

