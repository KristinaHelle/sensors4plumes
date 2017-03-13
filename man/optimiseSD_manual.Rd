\name{optimiseSD_manual}
\alias{optimiseSD_manual}
\title{
Plot cost map for interactive adding and deleting of sensors
}
\description{
Given a \code{costMap} function and a set of sensors it computes the cost map and plots it. Then the user can interactively add or delete sensors and plot the resulting cost and cost map..   
}
\usage{
optimiseSD_manual(simulations, costFun, costMap,
  locationsAll = 1:nLocations(simulations), locationsFix = integer(0), locationsInitial = integer(0),
  aimCost = NA, aimNumber = NULL,
  nameSave = NA, verbatim = FALSE,
  maxIterations = 10,
  valuesPlot = integer(0), colors = grey.colors)

}

\arguments{
Some arguments are the same for all optimisation algorithms, they are marked by a \code{*}, for detail see \code{\link{optimiseSD}}. Others are directly forwarded to \code{\link[genalg]{rbga.bin}}, they are marked by \code{**}
  \item{simulations}{\code{*}}
  \item{costFun}{\code{*}}
  \item{costMap}{A function to compute a cost value for each location,
  must return a list with \code{"costLocations"} (vector of length equal to \code{nLocations(simulations)}); needs to have the same parameters (\code{simulations}, \code{locations}) as \code{costFun} and may also be prepared by \code{replaceDefault} with \code{type = "costFun.optimiseSD"}.}
  \item{locationsAll}{\code{*}}
  \item{locationsFix}{\code{*}}
  \item{locationsInitial}{\code{*}}
  \item{aimCost}{\code{*} ignored}
  \item{aimNumber}{\code{*} ignored}
  \item{nameSave}{\code{*} ignored}
  \item{verbatim}{ ignored}
  \item{maxIterations}{maximal number of iterations to add or delete sensors}
  \item{valuesPlot}{names of values in the data associated with the locations of \code{simulations} to be plotted in addition to the cost map}
  \item{colors}{color ramp for plotting - a function like \code{grey.colors})}
}

%\details{
%
%}

\value{
A list
\item{cost}{cost in each iteration, based on fix and current sensors}
\item{identify}{Identity of points chosen in each iteration; these identities are 'raw' -for "SpatialIndexDataFrame" and "SpatialPolygridDataFrame" they may differ from the indices of the acutal locations as data is transformed into SpatialPixels before plotting, they may also include invalid and multiple choices.}
\item{SD}{indices of sampling design in each iteration; \code{locationsFix} are not included}
}
\author{ kristina.helle@uni-muenster.de  }

\examples{
# example 1: equal spread on "SpatialPolygrid"
# prepare data and functions
data(radioactivePlumes_local)
radioactivePlumes_loc = subset(radioactivePlumes_local, locations = 1:2001)

meanFun = function(x){mean(x, na.rm = TRUE)}
spatialSpreadMinDist = replaceDefault(
  spatialSpread,
  newDefaults = list(
    weightByArea = TRUE, 
    fun = minimalDistance, 
    fun_R = meanFun),
  type = "costFun.optimiseSD"
  )[[1]]

radioactivePlumes_loc@locations@data$p1 = getValues(
  subset.Simulations(radioactivePlumes_loc, plumes = 1, values = 1)@values)

# interactive optimisation
optSD_manual1 = optimiseSD_manual(simulations = radioactivePlumes_loc, 
                                  costFun = spatialSpreadMinDist, 
                                  costMap = spatialSpreadMinDist,
                                  locationsFix = seq(1, 2001, 100))


# example 2: detection on SpatialPolygons (small artificial data)
# prepare data (toy example) and functions
data(SPolygonsDF)
simulationsPolygons = subset(radioactivePlumes_local, 
                           locations = c(200, 221, 242), plumes = 31:35)
simulationsPolygons@locations = SPolygonsDF

simulationsPolygons@values$detectable = calc(
  simulationsPolygons@values$maxdose,
  fun = function(x){x >= 1e-3}) 
simulationsPolygons@plumes$totalDose = 
  summaryPlumes(simulationsPolygons, fun = sum, values = "finaldose")[[2]]

optSD_manual2 = optimiseSD_manual(simulations = simulationsPolygons, 
                                  costFun = singleDetection, 
                                  costMap = singleDetection,
                                  col = heat.colors)
}
