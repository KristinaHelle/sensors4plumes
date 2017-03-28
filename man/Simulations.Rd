\name{Simulations-class}
\docType{class}

\alias{Simulations}
\alias{Simulations-class}
\alias{coerce,Simulations,SpatialDataFrame-method}
\alias{nLocations}
\alias{nLocations.Simulations}
\alias{nLocations,Simulations-method}
\alias{nPlumes}
\alias{nPlumes.Simulations}
\alias{nPlumes,Simulations-method}
\alias{nKinds}
\alias{nKinds.Simulations}
\alias{nKinds,Simulations-method}


\title{Class "Simulations" }
\description{  Class for maps of plume simulations with spatial reference and plume properties. 
}

\section{Objects from the Class}{
Objects hold plume simulations and additional data: plume simulations are maps, all defined by values at the same set of locations; there may be several maps for each plume (e.g. the concentration of different pollutants in the same accident scenario, values at different times, etc.). 
As all maps have the same spatial properties, it is sufficient to hold the \code{locations} once as a \code{\link{SpatialDataFrame}}. 
The simulations themselves can be represented by a matrix where rows refer to location and columns refer to different plumes; in case of several maps per plume this can be extended to an array. These \code{values} arrays are actually implemented as \code{\link[raster]{raster}} objects; thus they can hold more data than it fits into memory. 
To keep information that refers to the single plumes like the kind of scenario they belong to or their probability to occur, a \code{\link{data.frame}} of \code{plumes} is part of each \code{Simulations} object.

%Objects consist of maps of several plumes (\code{values}), a \code{data.frame} with global properties of each plume (\code{plumes}) and a \code{SpatialDataFrame}  with properties of all locations where plume values are given (\code{locations}).
%The \code{values} are given as \code{raster} objects where each row belongs to a location and each column belongs to a plume - note that this means that raster layers are not maps. 
Objects may be created by calls to the function \code{Simulations(locations, plumes, values)}. 
}

\section{Slots}{
  \describe{
    \item{\code{locations}:}{\code{\link{SpatialDataFrame}}; locations and location-related data, e.g. population or cost to put a sensor here (n locations)}
    \item{\code{plumes}:}{\code{data.frame}; data related to plumes but not to locations, e.g. average impact of the plume (aggregated over space) or likelihood of this plume to occur  --- rows refer to plumes, columns to different kind of information (N plumes)}
    \item{\code{values}:}{\code{\link[raster]{raster}}; maps of plume values; each column represents the map of one plume: raster layers do not represent maps, but contain the values of the maps of all plumes (n x N values) --- spatial properties of the \code{values} are meanigless and must have the standardised form ---; there can be different layers to contain different kind of data, e.g. cumulative concentrations of a pollutant and time when pollution started.}
  }
}

%\section{Extends}{}

\section{Methods}{
\describe{
\item{nLocations}{ number of \code{locations} (equals \code{nrows(values)})
}
\item{nPlumes}{ number of \code{plumes} (equals \code{ncols(values)})
}
\item{nKinds}{ number of kind of \code{values} (equals \code{nlayers(values)})
}
%\item{subset}{ subsetting of \code{Simulations} by locations or plumes, there is no \code{"["} method; for details see \code{\link{subset.Simulations}}}
\item{cbind}{ cbind-like method to combine the plumes of different \code{Simulations} with same locations and value types; for details see \code{\link{cbind}}
}
\item{extractSpatialDataFrame}{ extract one \code{plume} from one layer of the \code{values} and assign it to the \code{locations} and return the resulting \code{SpatialDataFrame}; ; for details see \code{\link{extractSpatialDataFrame}}
}
}
}

\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

%\seealso{
%\code{\link{}}}

\examples{
# generate Simulations object: small, artificial example
data(SPixelsDF)
plumes = data.frame(source = c("A", "A", "B", "B", "B"),
                    date = c("2000-01-01", "2000-04-01", 
                             "2000-07-01", "2000-01-01", "2000-01-03"),
                    totalCost = runif(5, min = 5, max = 15))

values1 = replicate(n = nrow(plumes), expr = rlnorm(length(SPixelsDF), sdlog = 2))
values2 = replicate(n = nrow(plumes), expr = rnorm(length(SPixelsDF), m = 10, sd = 3))
values = stack(raster(x = values1, xmn = -90, xmx = 90, ymn = -90, ymx = 90, 
                      crs = "+init=epsg:4326" ), 
               raster(x = values2, xmn = -90, xmx = 90, ymn = -90, ymx = 90, 
                      crs = "+init=epsg:4326" ))

Simulations1 = Simulations(locations = SPixelsDF, plumes, values)

# nLocations, nPlumes, nKinds
nLocations(Simulations1)
nPlumes(Simulations1)
nKinds(Simulations1)

}
\keyword{classes}
