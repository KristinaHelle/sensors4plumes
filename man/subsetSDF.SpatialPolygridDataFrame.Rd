\name{subsetSDF.SpatialPolygridDataFrame}

\alias{subsetSDF.SpatialPolygridDataFrame}
\alias{subsetSDF,SpatialPolygridDataFrame-method}
%\alias{subSDF.SpatialPolygridDataFrame}

\title{ Subsetting objects of class SpatialPolygridDataFrame}

\description{ Subsetting of the grid by either indices or coordinates, or by a logical mask; subsetting the index or the data. }

\usage{
subsetSDF.SpatialPolygridDataFrame(x, locations, data = names(x@data), ..., 
  coord_x, coord_y, grid_i, grid_j, grid_ij)                
}

\arguments{
\item{x}{ object of class \link{SpatialPolygridDataFrame} whereof we want a subset}
\item{locations}{ integer vector: the locations of x with these indices to be kept (refers to \code{x@index}); order is taken into account, multiplicity is ignored }
\item{data}{ integer vector or names of columns of the data of x to be kept (refers to \code{x@data}) }
\item{...}{further arguments to be passed to other methods}
\item{grid_i}{ integer vector: the columns of the grid of x to be kept, from top (North) to bottom (South)  }
\item{grid_j}{ integer vector: the rows of the grid of x to be kept, from left (West) to right (East) }
\item{coord_x}{ two numerical values: min and max coordinates in x-direction (columns, from West), keep the part in between, including the grid cells with these x-coordinates}
\item{coord_y }{two numerical values: min and max coordinates in y-direction (rows, from South), keep the part in between, including the grid cells with these y-coordinates }
\item{grid_ij}{matrix of logical, same size as x (ncol(grid_ij) must be \code{x@grid@grid.cells.dim[1]}, \code{nrow(grid_ij)} must be \code{x@grid@grid.cells.dim[2]}); keep values in grid cells where \code{grid_ij} is  \code{TRUE}}
For \code{grid_i}, \code{grid_j} order is ignored, original order of grid is kept.
}

\value{
%A \emph{SpatialPolygridDataFrame}; the subsetting can be done cumulatively, e.g. if \code{index_i} and \code{coord_i} is given, only the colums that have an index of \code{index_i} and coordinates within the limits of \code{coord_i} are returned. The result is always a complete, regular grid, nonselected parts within the grid are replaced by NA in the index. If, for example, \code{x} has 4 columns, then \code{subset(x, grid_i = c(1,3))} has 3 columns: columns 1:3 of \code{x}, column 2 is replaced by NA. It is determined separately for \code{grid_i} and \code{grid_j}, \code{coord_i} and \code{coord_j}, \code{grid_ij}, and \code{grid_ij} how much the grid is cropped, the minimal common grid is chosen -- see examples below.
All subsetting parameters can be combined, the result is the intersection of all parameters, e.g. if \code{index_i} and \code{coord_y} is given, only the colums that have an index of \code{index_i} and coordinates within the limits of \code{coord_y} are returned. 
If invalid values are chosen (outside of the grid), a warning is thrown and the intersection of the parameters with \code{x} is used. -- Only for \code{data}, invalid values result in an error.

The result is a SpatialPolygridDataFrame, cropped to the extent of the remaining grid cells with values. Nonselected parts within this grid are replaced by NA in the index. 
%If the resulting grid has only one row or column, the result is provided as SpatialPointsDataFrame with points at the centers of all remaining grid cells and the associated values as data.

%If, for example, \code{x} has 4 columns, then \code{subset(x, grid_i = c(1,3))} has 3 columns: columns 1:3 of \code{x}, column 2 is replaced by NA. It is determined separately for \code{grid_i} and \code{grid_j}, \code{coord_i} and \code{coord_j}, \code{grid_ij}, and \code{grid_ij} how much the grid is cropped, the minimal common grid is chosen -- see examples below.. 
}

%\note{}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }
%\seealso{ \link{}}

\examples{
data(SPolygridDF)
# subset only attributes
SPolygridDF1 = subsetSDF(SPolygridDF, data = 2)
# subset by locations
SPolygridDF2 = subsetSDF(SPolygridDF, locations = c(1,3,5,7))
# subset by rows and columns of the grid
SPolygridDF3 = subsetSDF(SPolygridDF, grid_i = 3:4, grid_j = 3:9)
# subset by individual grid cells
Grid_ij = matrix(nrow = 6, ncol = 6, byrow = TRUE, data = c(
                         TRUE, FALSE, FALSE, FALSE, FALSE, TRUE,
                        FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                        FALSE, FALSE,  TRUE,  TRUE, FALSE, FALSE,
                        FALSE, FALSE,  TRUE,  TRUE, FALSE, FALSE,
                         TRUE, FALSE, FALSE, FALSE, FALSE, TRUE,
                        FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
SPolygridDF4 = subsetSDF(SPolygridDF, grid_ij = Grid_ij)
# subset by bounding box of coordinates
SPolygridDF5 = subsetSDF(SPolygridDF, coord_x = c(3,11), coord_y = c(3,11))

# combined subset (overlap of all subsetting parameters), result has 1 cell
SPolygridDF6 = subsetSDF(SPolygridDF, 
data = "a", locations = 4:1, coord_x = c(0,6), grid_i = 1:3)
}

%\keyword{}
