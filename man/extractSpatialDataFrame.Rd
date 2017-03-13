\name{extractSpatialDataFrame}

\alias{extractSpatialDataFrame}
%\alias{extractSpatialDataFrame.Simulations}
%\alias{extractSpatialDataFrame,Simulations-method}

\title{
Extract some values of Simulations to a SpatialDataFrame with the same spatial properties.
}
\description{
Turn some of the \code{values} of a \code{Simulations} object into a \code{\link{SpatialDataFrame}}, keeping the spatial reference.
}
\usage{
extractSpatialDataFrame(obj, kinds = 1:nlayers(obj@values),
  plumes = 1:ncol(obj@values), chunksize = 1e+7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{\code{\link{Simulations}} object
  }
  \item{kinds}{integer: values of which layers (kind) of the \code{values} to extract 
  }
  \item{plumes}{integer: values of which plume to extract 
  }
  \item{chunksize}{limit of cells to extract, if selected kinds and plumes exceed this limit, only the first are extracted (with warning)}
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
