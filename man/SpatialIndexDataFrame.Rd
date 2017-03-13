\name{SpatialIndexDataFrame-class}
\docType{class}

\alias{SpatialIndexDataFrame}
\alias{SpatialIndexDataFrame-class}
\alias{SpatialIndexDataFrame-method}
\alias{coerce,data.frame,SpatialIndexDataFrame-method}
\alias{coordinates,SpatialIndexDataFrame-method}
\alias{proj4string,SpatialIndexDataFrame-method}
\alias{bbox,SpatialIndexDataFrame-method}
\alias{spplot,SpatialIndexDataFrame-method}
\alias{length,SpatialIndexDataFrame-method}
\alias{cbind.SpatialIndexDataFrame}
\alias{rbind.SpatialIndexDataFrame}


\title{Class "SpatialIndexDataFrame" }

\description{  
This class is an extension of the \link[base]{data.frame} class. In addition it has an index vector of arbitrary length, consisting of integers that assign the rows of the data to the positions in the index.
The class is needed to complete the \link{SpatialDataFrame-class} in cases where no spatial information is available.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{as(x,"SpatialIndexDataFrame")}, from \code{x} of class \code{data.fame}, or from scratch with \code{SpatialIndexDataFrame(index, data.frame)}. 
}

\section{Slots}{
\describe{
\item{\code{data}:}{
Object of class \code{data.frame}.
}
\item{\code{index}:}{
Object of class \code{integer}; it must consist of the values \code{1:nrow(data)} and may contain \code{NA}. The length is arbitrary.
}
\item{\code{bbox}:}{
Object of class \code{matrix}; no values (only needed to fit the structure of \code{SpatialDataFrame}).
}
\item{\code{proj4string}:}{
Object of class \code{\link[sp]{CRS}}; no values (only needed to fit the structure of \code{SpatialDataFrame}).
}
}
}

%\section{Extends}{}

\section{Methods}{
\describe{
\item{coordinates}{
Returns \code{NULL}, only needed to fit the structure of the \code{SpatialDataFrame} class.
}
\item{proj4string}{
Returns \code{NA}, only needed to fit the structure ofthe \code{SpatialDataFrame} class.
}
\item{bbox}{
Returns  matrix of 0s, only needed to fit the structure ofthe \code{SpatialDataFrame} class.
}
\item{spplot}{
Turns index into a column (top to bottom) of grid cells with the respective values and plots it by \code{\link[sp]{spplot}}. Layout parameters like 
colours may be forwarded to spplot.
}
\item{cbind}{
Method to combine values of \code{SpatialIndexDataFrame} objects with identical \code{index}.
}
\item{rbind}{
Method to combine \code{SpatialIndexDataFrame} objects with \code{data} of identical type and name.
}
\item{subsetSDF}{See \link{subsetSDF.SpatialIndexDataFrame}. There is no \code{"["} function for this class. 
}
\item{areaSDF}{
See \link{areaSDF}, returns vector of 0s, needed to fit the structure of the \code{SpatialIndexDataFrame} class.
}
\item{length}{
\code{signature(x = "SpatialIndexDataFrame")}:
number of elements, i.e. of groups of cells with respectively share the attributes.
}
}
}

\seealso{
  \code{\link{SpatialDataFrame-class}}, which is a wrapper.
}

\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }


\examples{
# generate SpatialIndexDataFrame from scratch
index = c(2,2,1,1,1,1,3,2,3,3,3,1)
SIndexDF1 = SpatialIndexDataFrame(index = as.integer(index), 
                          data = data.frame(a = c(1,2,4), 
                                            b = c(0.1, 0.2, 0.3), 
                                            c = c("A", "B", "A")))
# generate from data.frame
data(USArrests)
SIndexDF2 = as(USArrests,"SpatialIndexDataFrame")


# spplot
spplot(SIndexDF1, zcol = c("b", "c"), col.regions = grey.colors(10))

# cbind
SIndexDF3 = cbind(SIndexDF1, SIndexDF1)
}
%\keyword{}
