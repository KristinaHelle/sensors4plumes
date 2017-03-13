\name{spatialSpreadFunctions}
\alias{minimalDistance}
\alias{krigingVariance}

\title{
Cost functions dependent only on sensor locations
}
\description{
\code{Kriging variance} computes the ordinary kriging variance at a given set of points for the case of certain sensor locations and variogram model.

\code{minimalDistance} computes for each given point the distance to the next sensor.

Both functions can be used as input to \code{spatialSpread} that runs them and turns the result into global cost.
}
\usage{
  krigingVariance(allLocations, locations, model)
  minimalDistance(allLocations, locations, algorithm = "kd_tree")
}

\arguments{
  \item{allLocations}{ \code{SpatialDataFrame} (not "SpatialIndexDataFrame"), for each of the locations a cost is computed}
  \item{locations}{indices of sensors in \code{allLocations}}
  \item{model}{\code{model} to be forwarded to \code{\link[gstat]{krige}}
  }
  \item{algorithm}{\code{character} to be forwarded as \code{algorithm} to \code{\link[FNN]{get.knnx}}}
}

\value{
A vector of length \code{allLocations} with cost for these locations.
}

\author{Kristina B. Helle, \email{kristina.helle@uni-muenster.de}}

\examples{
# see spatialSpread
}

