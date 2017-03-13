\name{subsetSDF.SpatialIndexDataFrame}

\alias{subsetSDF.SpatialIndexDataFrame}
\alias{subsetSDF,SpatialIndexDataFrame-method}

\title{ Subsetting objects of class SpatialIndexDataFrame}

\description{ Subsetting of objects of class \code{SpatialIndexDataFrame}: selects index and/or attribute(s), and returns an object of class \code{SpatialIndexDataFrame}. }

\usage{             
subsetSDF.SpatialIndexDataFrame(x, locations, data = names(x@data), ..., grid)
}
 
\arguments{
\item{x}{ object of class \code{SpatialIndexDataFrame}}
\item{locations}{ integer vector: keep data where \code{x@index} has this value (multiplicity ignored, order taken into account) }
\item{data}{ integer vector or names of columns of the \code{x@data} to be kept (order and multiplicity ignored)}
\item{...}{further arguments to be passed to other methods}
\item{grid}{ integer vector: indices of \code{x@index}, keep the corresponding data }
}


\value{
A \code{SpatialIndexDataFrame}; the subsetting can be done cumulatively, e.g. if \code{grid} and \code{locations} is given, only the index entries that fullfill both criteria are kept.
}

%\note{}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }
%\seealso{ \link{}}

\examples{
data(SIndexDF)
# subset
SIndexDF1 = subsetSDF(SIndexDF, grid = which(SIndexDF@index == 1)) 
SIndexDF2 = subsetSDF(SIndexDF, locations = 1) # identical to x_1
SIndexDF3 = subsetSDF(SIndexDF, grid = c(2:4, 6:8, 10:12, 14:16, 18:20), data = "c")
}

%\keyword{}
