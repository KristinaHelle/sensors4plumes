\name{SpatialDataFrame-class}

\docType{class}

\alias{SpatialDataFrame}
\alias{SpatialDataFrame-class}
%\alias{SpatialPolygridDataFrame}
%\alias{SpatialIndexDataFrame}
%\alias{SpatialPointsDataFrame}
%\alias{SpatialPixelsDataFrame}
%\alias{SpatialPolygonsDataFrame}

%\alias{subsetSDF}
%\alias{areaSDF}
\alias{coordinates,SpatialDataFrame-method}
\alias{proj4string,SpatialDataFrame-method}
\alias{bbox,SpatialDataFrame-method}
\alias{spplot,SpatialDataFrame-method}
\alias{length,SpatialDataFrame-method}
\alias{is.SpatialDataFrame}

\title{Class "SpatialDataFrame" }
\description{  Class for a data frame with associated spatial attributes (points, lines, polygons, (irregular) grid). }

\usage{
is.SpatialDataFrame(SpatialDataFrame)
}
\arguments{
\item{SpatialDataFrame}{an object, if its class belongs to one of the subclasses of \code{SpatialDataFrame-class} the function returns \code{TRUE} else \code{FALSE}}
}

\section{Objects from the Class}{
This class is a wrapper for some classes from package \code{\link[sp]{sp}}: 

\link[sp]{SpatialPointsDataFrame-class},

\link[sp]{SpatialPixelsDataFrame-class},

%\link[sp]{SpatialGridDataFrame-class},

\link[sp]{SpatialPolygonsDataFrame-class},

\link[sp]{SpatialLinesDataFrame-class}

Besides it contains two classes specific to \code{sensors4plumes}: 

\link{SpatialIndexDataFrame-class}: Objects of this class have no spatial reference, they are just \code{data.frames} with an extra index.

\link{SpatialPolygridDataFrame}: These objects are similar to those from \link[sp]{SpatialGridDataFrame-class}, but several grid cells may refer to the same attributes in the data frame.
}



%\section{Extends}{}


\section{Methods}{

The following methods exist for all objects of class \emph{SpatialDataFrame}. They may be defined individually for the respective subclasses, refer to the subclasses for details.
\describe{
\item{is.SpatialDataFrame}{\code{TRUE} for objects of class \code{SpatialDataFrame}, else \code{FALSE}}
\item{length}{\code{signature(x = "SpatialDataFrame")}: length (number of spatial elements)
}
\item{coordinates}{\code{signature(x = "SpatialDataFrame")}: retrieves coordinates (of centroids)
}
\item{proj4string}{\code{signature(x = "SpatialDataFrame")}: retrieves the projection as \link[sp]{CRS-class}
}
\item{bbox}{\code{signature(x = "SpatialDataFrame")}: retrieves the bounding box
}
\item{spplot}{\code{signature(x = "SpatialDataFrame", zcol = "integer", ...)}: plots maps of the attributes \code{zcol} ("integer" or names of \code{x@data}); wrapper of \code{\link[sp]{spplot}}
}
%\item{cbind}{\code{signature(...)}: binds columns of the data if all other slots are identical, else it throws an error.}
\item{subsetSDF}{\code{signature(x = "SpatialDataFrame", locations = "integer", data = "integer")}: selects spatial sub-objects (\code{locations}) and/or attributes (\code{data}: "integer" or "character" names of \code{x@data}), returns an  object of the same class. This is a general method for all \code{SpatialDataFrame} objects. There are special \code{subset} or \code{"["} functions for the subclasses that may have more possibilities to choose (e.g. subset grid by a bounding box), see the subclasses for details.
}        
\item{areaSDF}{\code{signature(x = "SpatialDataFrame")}: computes the area/length of the elements (grid cells, polygons...), for point data it returns 0s.
}
}
Refer to the subclasses for further methods.  
}


\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

%\seealso{
%For examples also refer to the subclasses.
%}

\examples{
# generate a SpatialPixelsDataFrame
data(meuse.grid)
gridded(meuse.grid) = ~ x + y

# extract some locations and attibutes
meuse.subgrid = subsetSDF(meuse.grid, 
        locations = c(1:45), 
        data = c("dist", "soil", "ffreq"))

# plot
spplot(meuse.subgrid, zcol = 1:2)
}
\keyword{classes}
