\name{SpatialPolygridDataFrame-class}
\docType{class}

\alias{SpatialPolygridDataFrame}
\alias{SpatialPolygridDataFrame-class}
\alias{SpatialPolygridDataFrame-method}
\alias{coerce,SpatialGridDataFrame, SpatialPolygridDataFrame-method}
\alias{coerce,SpatialPolygridDataFrame, SpatialGridDataFrame-method}
\alias{coerce,SpatialPointsDataFrame, SpatialPolygridDataFrame-method}
\alias{coordinates,SpatialPolygridDataFrame-method}
\alias{proj4string,SpatialPolygridDataFrame-method}
\alias{bbox,SpatialPolygridDataFrame-method}
\alias{spplot,SpatialPolygridDataFrame-method}
\alias{length,SpatialPolygridDataFrame-method}
\alias{cbind.SpatialPolygridDataFrame}


\title{Class "SpatialPolygridDataFrame"}

\description{Class for a data frame with associated spatial grid, several grid cells may share one attribute.}

\section{Objects from the Class}{
Objects can be created by calls to the function \code{SpatialPolygridDataFrame}.
They may also be created from \code{\link[sp]{SpatialGridDataFrame}} objects by calls of the form \cr \code{y = as(x, "SpatialPolygridDataFrame")}, in this case all grid cells have individual attributes. 
Transformation in the opposite direction can be done by \code{x = as(y, "SpatialGridDataFrame")}; it copies the attributes to all assigned grid cells.
SpatialPolygridDataFrame objects can also be created from \code{\link[sp]{SpatialPointsDataFrame}} objects by calls of the form \code{as(x, "SpatialPolygridDataFrame")};  it generates a rectangular grid that covers all points and then assigns cells to the nearest neighbouring point -- it may fail if points are too irregular to allow for a common grid.
}

\section{Slots}{
\describe{
\item{\code{data}:}{
object of class \code{\link[base]{data.frame}}; contains the data, one row per index 
}
\item{\code{index}:}{
object of class \code{\link[base]{integer}}, assigns data to grid cells:
order indicates reference to grid cells: upper left to right, then to bottom -- length must agree with number of grid cells; 
value indicates the row of the \code{data} with the related attributes; may contain \code{NA}, and values outside of \code{1:nrow(data)} are ignored.
}
\item{\code{grid}:}{
object of class \code{\link[sp]{GridTopology}}
}
\item{\code{bbox}:}{
bounding box of the grid
}
\item{\code{proj4string}}{
object of class \code{\link[sp]{CRS}}; projection of the grid, may be missing
}
}
}

\section{Methods}{
\describe{
\item{coordinates}{
\code{signature(obj = "SpatialPolygridDataFrame")}: 
Retrieves the mean coordinates for each index, i.e. the centroid of all grid cells belonging to it. 
}
\item{proj4string}{
\code{signature(obj = "SpatialPolygridDataFrame")}: 
Returns the projection character string, the value of the slot \code{proj4string}.
}
\item{bbox}{
\code{signature(obj = "SpatialPolygridDataFrame")}: 
Returns the bounding box matrix, the value of the slot \code{bbox}.
}
\item{spplot}{
\code{signature(obj = "SpatialPolygridDataFrame", geoTiffPath, zcol = names(obj@data), plot = TRUE, returnSGDF = FALSE ...)}: 
Method to plot maps, by coercing the \code{SpatialPolygridDataFrame} into a \code{\link[sp]{SpatialGridDataFrame}} -- see above -- and using \code{\link[sp]{spplot}} to plot it. Therefore it may take any layout parameters to be forwarded to \code{spplot}.
To keep the resulting \code{SpatialGridDataFrame} it may be returned by \code{returnSGDF = TRUE}. It can also be saved to a (multilayer) tif file by providing a file path in \code{geoTiffPath}. 
To speed up the function, it may be run without graphical output by \code{plot = FALSE}, e.g. to directly create a tif file.
}
\item{cbind}{
\code{signature(...)}: 
\code{cbind}-like method; if grids and indices agree, the data are combined by cbind. 
}
\item{subsetSDF}{
Subsetting, see \code{\link{subsetSDF,SpatialPolygridDataFrame-method}}; there is no \code{"["} function.
}
\item{areaSDF}{
Method to compute the areas of all polygrids, see \code{\link{areaSDF,SpatialPolygridDataFrame-method}}.
}
\item{length}{
\code{signature(x = "SpatialPolygridDataFrame")}:
number of elements, i.e. of groups of cells with respectively share the attributes.
}
}
}

\note{
The code for coercing is based on \code{\link{polygrid2grid}} and \code{\link{points2polygrid}}. Using these functions directly allows better control of the parameters, like the resolution of the grid.
}

\seealso{
  \code{\link{SpatialDataFrame-class}}, which is a wrapper.
}

\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
# prepare 
index1 = as.integer(
  c( 6, 6, 7, 7, 8, 8,
     6, 6, 7, 7, 8, 8,
     5, 5, 1, 2, 9, 9,
     5, 5, 4, 3, 9, 9,
     12,12,11,11,10,10,
     12,12,11,11,10,10))
dataFrame1 = data.frame(a = 1:12 * 10, b = 2^(12:1))
grid1 = GridTopology(c(1,1), c(2,2), c(6,6))
spatialGrid1 = SpatialGrid(grid1, CRS("+proj=longlat +datum=WGS84"))

#- - create object by function SpatialPolygridDataFrame - -#
SPolygridDF1 = SpatialPolygridDataFrame(
  grid = grid1,
  data = dataFrame1,
  index = index1,
  proj4string = CRS("+proj=longlat +datum=WGS84"))

#- - - - coerce from SpatialGridDataFrame - - - - - - - -#
spatialGridDataFrame1 = SpatialGridDataFrame(
    grid = spatialGrid1, data = dataFrame1[index1,])
# coerce, each cell keeps individual values
s = as(spatialGridDataFrame1, "SpatialPolygridDataFrame")

# coerce back (name of coordinates may change)
#spatialGridDataFrame2 = as(SPolygridDF2, "SpatialGridDataFrame")
spatialGridDataFrame2 = as(s, "SpatialGridDataFrame")

#- - - - coerce from SpatialPointsDataFrame - - - - - - - #
# from irregular points, grid is created, cells are assigned values of nearest point
spatialPoints1 = SpatialPoints(coordinates(SPolygridDF1) + runif(24))
spatialPointsDataFrame1 = SpatialPointsDataFrame(
    coords = spatialPoints1, data = dataFrame1)
\dontrun{
## takes some seconds
SPolygridDF3 = as(spatialPointsDataFrame1, "SpatialPolygridDataFrame")


# from regular points: one cell per point
spatialPointsDataFrame2 = as(spatialGridDataFrame1, "SpatialPointsDataFrame")
SPolygridDF4 = as(spatialPointsDataFrame2, "SpatialPolygridDataFrame")

# - - - - - - - - - - - - - coordinates - - - - - - - - - - - - #
# compare irregular input points and centroids of output
plot(spatialPoints1)
points(coordinates(SPolygridDF3))

# - - - - - - - - - - - - - proj4string - - - - - - - - - - - -  #
# only to retrieve projection attributes
#proj4string(SPolygridDF2)
proj4string(s)
proj4string(SPolygridDF3)
# - - - - - - - - - - - - - proj4string - - - - - - - - - - - -  #
bbox(SPolygridDF3)
# - - - - - - - - - - - spplot - - - - - - - - - - - - - - - - #

## takes some time
# plot map with original points
spplot(SPolygridDF3, zcol = "a", 
  sp.layout = list(list("sp.points", spatialPoints1, col = 8)))
# plot map with centres of polygrid cells
spplot(SPolygridDF4, zcol = "a", col.regions = rainbow(16),
  sp.layout = list("sp.points", spatialPointsDataFrame2@coords, col = 8))


# - - - - - - - - - - - - cbind - - - - - - - - - - - - - - - - #
#SPolygridDF5 = cbind(SPolygridDF2, SPolygridDF4) 
SPolygridDF5 = cbind(s, SPolygridDF4) 
}
}
\keyword{classes}
