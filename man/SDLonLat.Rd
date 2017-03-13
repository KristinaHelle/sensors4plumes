\name{SDLonLat}
\alias{SDLonLat}
\title{Transform locations of simulations into longitude latitude coordinates
}
\description{
Given a \code{\link{Simulations}} object and indices of locations it returns \code{\link[sp]{SpatialPoints}} in longitude and latitude.}

\usage{
SDLonLat(simulations, SD)
}
\arguments{
  \item{Simulations}{a \code{\link{Simulations}} object}
  \item{SD}{vector of locations indices or matrix (with SDs in rows) or list}
}
\value{
\code{List of \link[sp]{SpatialPoints} with WGS84 longlat projection.}
}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
data(radioactivePlumes_area)
SDs1 = matrix(sample.int(2500, 40), nrow = 8)

LL1 = SDLonLat(simulations = radioactivePlumes_area, 
                   SD = SDs1)
}

