\name{spatialSpread}
\alias{spatialSpread}

\title{
Cost function based on sensor locations only
}
\description{
This is a wrapper function to compute cost of a set of sensors based on their coordinates only. It takes as parameter \code{fun} an actual function to compute such cost for each location of the provided \code{simulations} and can summarise it with a function \code{fun_R}
}
\usage{
spatialSpread(simulations, locations, weightByArea = TRUE, fun, fun_R)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simulations}{ \code{Simulations} object or \code{SpatialDataFrame} to compute local cost in each \code{location} (\code{values} are ignored)}
  \item{locations}{ indices of \code{locations} of sensors }
  \item{weightByArea}{ if cost in the locations has to be multiplied by the area associated to them before summary; for \code{"SpatialPointsDataFrame"} the factor is 0.}
  \item{fun}{\code{function} to compute location-wise cost, for examples see \code{spatialSpreadFunctions}; needs the parameters \code{allLocations}, \code{locations}; can be the result of \code{replaceDefault} with \code{type = "fun.spatialSpread"}; result must have length number of locations}
  \item{fun_R}{\code{function} to summarise the result of \code{fun}; must have parameter \code{x}, generate by \code{replaceDefault} with \code{type = "funR.spatialSpread"}}
}

\details{
The function itself is only a wrapper to turn \code{simulations} and \code{locations} into the correct form and apply the functions.
}

\value{
List
 \item{cost}{result of \code{fun_R}, missing if no \code{fun_R} given}
 \item{costLocations}{result of \code{fun}, length equals number of locations in \code{simulations}}
}

\author{Kristina B. Helle, \email{kristina.helle@uni-muenster.de}}


\examples{
# prepare data and functions
data(radioactivePlumes)
data(medianVariogram)

krigingVarianceMedian = 
  replaceDefault(krigingVariance, newDefaults = list(model = medianVariogram))[["fun"]]

meanFun = function(x){mean(x, na.rm = TRUE)}

locationsSensors = sample.int(nLocations(radioactivePlumes), 50)

spatialSpread_minDist = spatialSpread(
    simulations = radioactivePlumes,
    locations = locationsSensors, 
    weightByArea = TRUE,
    fun = minimalDistance,
    fun_R = meanFun
)
spatialSpread_krigingVar = spatialSpread(
    simulations = radioactivePlumes,
    locations = locationsSensors, 
    weightByArea = TRUE,
    fun = krigingVarianceMedian,
    fun_R = meanFun
  )

# plot maps
\dontrun{
## takes some seconds
y = radioactivePlumes@locations
y@data$minDist = spatialSpread_minDist[["costLocations"]]
y@data$krigVar = spatialSpread_krigingVar[["costLocations"]]
yPoints = as(y, "SpatialPointsDataFrame")

# distance to next sensor
spplot(y, zcol = "minDist", 
       sp.layout = list("sp.points", yPoints[locationsSensors,], 
                          col = 3))
# kriging variance
spplot(y, zcol = "krigVar", 
       sp.layout = list("sp.points", yPoints[locationsSensors,], 
                          col = 3))
}
}
