\name{SDF2simulations}

%\alias{SDF2simulations.SpatialDataFrame}
\alias{SDF2simulations}
%\alias{SDF2simulations,SpatialDataFrame-method}

\title{
Turn a SpatialDataFrame into Simulations.
}
\description{
Turn a \code{\link{SpatialDataFrame}} into \code{Simulations}: turn some of the data into the \code{values}, keep the others as the data of the \code{locations}.
}
\usage{
SDF2simulations(x, indices = matrix(1:ncol(x@data), nrow = 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{\link{SpatialDataFrame}} 
  }
  \item{indices}{matrix of integers indicating the columns of the \code{data} to be transformed to \code{values}; each row stands for a \code{kind} of \code{values}
  }
}

\details{
This function can be used to extract the maps in \code{Simulations} objects: the values are assigned to the spatial properties, for plotting, combining it with other spatial data etc. As the size of all \code{values} may exceed the memory, it may extract only the first of the kinds and plumes. 
}
\value{
A \code{\link{SpatialDataFrame}} with the same spatial reference as the \code{locations} of the \code{obj} with the values of the chosen \code{plume} and \code{layer} as attribute.
}

\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
data(SimulationsSmall)
SimulationsSmall_134_21 = extractSpatialDataFrame(
  SimulationsSmall,
  plumes = c(1,3,4), kinds = c(2,1))
}
