\name{points2polygrid}

\alias{points2polygrid}

\title{ Turn points (and data) into a SpatialPolygridDataFrame}


\description{ \code{points2polygrid} determines a fine grid, such that each point is represented by a set of grid cells (nearest neighbourhood); to each grid cell it assigns the values of the nearest point. It may also assign cells of a given grid to the nearest neighbour points.}


\usage{
points2polygrid(points, grid, index, data, 
tolerance = signif(mean(diff(t(bbox(points))))/1000, digits = 1))           
}

\arguments{
\item{points}{ coordinates of the original points: \code{data.frame}, \code{matrix} with coordinates in columns 1 and 2, or \code{\link[sp]{SpatialPoints}} }
\item{grid}{ target grid: \code{\link[sp]{GridTopology}} or \code{\link[sp]{SpatialGrid}} }
\item{index}{ target index indicating for each grid cell the associated data: integer, length must fit number of cells in grid
(order: top left to right, then to bottom) }
\item{data}{ data associated with points, to be assigned to grid cells: \code{\link{data.frame}} or vector of same length as \code{points} }
\item{tolerance}{ lower limit of the resolution of the generated grid, default generates about 1000 x 1000 grid cells; if points lie on a regular grid, the resulting resolution is automatically adapted}
}

\details{
It is sufficient to provide \code{points}, they can be used to generate a grid and associate the cells to the closest points. The generated grid covers all points; resolution is chosen such that there are cell boundaries between the x-coordinates and y-coordinates of all points respectively.  If \code{points} lie on a regular grid, resolution is chosen to fit this grid.
A \code{grid} can also be provided which saves computational effort. 

If no \code{index} is given, grid cells are assigned to the nearest point (based on \code{\link[FNN]{get.knnx}}).
If \code{index}, \code{grid} and \code{points} are given, \code{points} are ignored.
If \code{points} are at the same location or so close that they are merged, the associated data of only one point is used. 
}

\value{
A \code{\link{SpatialPolygridDataFrame}}.
If no \code{data} is given, the returned data consists of one column 'Index' with the indices of the associated points.
\code{proj4string} of \code{points} is kept if available, else the one of \code{grid} is taken.
}

\note{
This function is used to coerce \emph{SpatialPointsDataFrame} to \emph{SpatialPolygridDataFrame}.
}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

%\seealso{ \link{}}

\examples{
# prepare
data(SPolygridDF)
regularPoints = coordinates(SPolygridDF)
irregularPoints = regularPoints + runif(12, 0.1)

# grid for regular points
newSPolygrid_reg = points2polygrid(points = regularPoints, 
                                   data = data.frame(a = 12:1 * 5))
spplot(newSPolygrid_reg, 
       sp.layout = list(list("sp.points", regularPoints, col = 3)))
       
# generate grid for irregular points with given resolution
newSPolygrid_irreg = points2polygrid(points = irregularPoints, tolerance = 0.5)
spplot(newSPolygrid_irreg, sp.layout = list("sp.points", irregularPoints, col = 3))

# with given grid
newSPolygrid_grid = points2polygrid(points = irregularPoints,
  grid = GridTopology(c(2,2), c(1,1), c(10,10)))
spplot(newSPolygrid_grid, sp.layout = list("sp.points", irregularPoints, col = 3))

}

%\keyword{}
