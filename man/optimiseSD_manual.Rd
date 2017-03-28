\name{optimiseSD_manual}
\alias{optimiseSD_manual}
\title{
Plot cost map for interactive adding and deleting of sensors
}
\description{
Given a \code{costMap} function and a set of sensors it computes the cost map and plots it. Then the user can interactively add or delete sensors and plot the resulting cost and cost map..
}
\usage{
optimiseSD_manual(simulations, costFun,
  locationsAll = 1:nLocations(simulations), locationsFix = integer(0),
  locationsInitial = integer(0),
  aimCost = NA, aimNumber = NA,
  nameSave = NA, plot = TRUE, verbatim = FALSE,
  costMap = NA, maxIterations = 10, 
  valuesPlot = integer(0), colors = grey.colors)
}

\arguments{
Some arguments are the same for all optimisation algorithms, they are marked by a \code{*}, for detail see \code{\link{optimiseSD}}. Others are directly forwarded to \code{\link[genalg]{rbga.bin}}, they are marked by \code{**}
  \item{simulations}{\code{*}}
  \item{costFun}{\code{*} ignored - cost is derived from \code{costMap}}
  \item{locationsAll}{\code{*}}
  \item{locationsFix}{\code{*}}
  \item{locationsInitial}{\code{*}}
  \item{aimCost}{\code{*} ignored}
  \item{aimNumber}{\code{*} ignored}
  \item{nameSave}{\code{*} ignored}
  \item{plot}{\code{*} ignored}
  \item{verbatim}{ ignored}
  \item{costMap}{A function to return a list of \cr
  \code{"cost"}: cost in general, a single value,  \cr
  \code{"costLocations"}: a cost value for each location(vector of length equal to \code{nLocations(simulations)}). \cr
  Needs to have the same parameters (\code{simulations}, \code{locations}) as \code{costFun} and may also be prepared by \code{replaceDefault} with \code{type = "costFun.optimiseSD"}.}

  \item{maxIterations}{maximal number of iterations to add or delete sensors}
  \item{valuesPlot}{names of values in the data associated with the locations of \code{simulations} to be plotted in addition to the cost map}
  \item{colors}{color ramp for plotting - a function like \code{grey.colors})}
}

%\details{
%
%}

\value{
A list, the first two entries are common to all optimisation algorithms, they are marked with \code{*}, see \code{\link{optimiseSD}} for details.
\item{SD}{\code{*} best sampling designs}
\item{evaluation}{\code{*} cost and size of \code{SD}}
\item{report}{ a list of \cr
\code{cost}:  cost in each iteration, based on fix and current sensors \cr
\code{identify}:  identity of points chosen in each iteration; these identities are 'raw' -for "SpatialIndexDataFrame" and "SpatialPolygridDataFrame" they may differ from the indices of the actual locations as data is transformed into SpatialPixels before plotting, they may also include invalid and multiple choices.\cr
\code{locationsCurrent}: indices of sampling design in each iteration; \code{locationsFix} are not included}
}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
# prepare data and functions
data(radioactivePlumes)

meanFun = function(x){mean(x, na.rm = TRUE)}
spatialSpreadMinDist = replaceDefault(
  spatialSpread,
  newDefaults = list(
    weightByArea = TRUE,
    fun = minimalDistance,
    fun_R = meanFun),
  type = "costFun.optimiseSD"
  )[[1]]

radioactivePlumes@locations@data$p1 = getValues(
  subset(radioactivePlumes, plumes = 1, kinds = 1)@values)

optimSD_man_minDist = replaceDefault(
  optimiseSD_manual,
  newDefaults = list(
    costMap = spatialSpreadMinDist
    )
)[["fun"]]

\dontrun{
## interactive optimisation
# inside optimiseSD
optSD_manual1 = optimiseSD(simulations = radioactivePlumes,
                          costFun = spatialSpreadMinDist,
                          optimisationFun = optimSD_man_minDist,
                          locationsFix = seq(1, 2001, 300),
                          locationsInitial = seq(1, 2001, 70),
                          locationsAll = setdiff(1:2001, seq(1,2001, 30)))

# directly using optimiseSD_manual
optSD_manual2 = optimiseSD_manual(simulations = radioactivePlumes,
                                  costFun = spatialSpreadMinDist,
                                  costMap = spatialSpreadMinDist,
                                  locationsFix = seq(1, 2001, 300),
                                  locationsInitial = seq(1, 2001, 70),
                                  locationsAll = setdiff(1:2001, seq(1,2001, 30)))
}
## the result of such a manual optimisation is in data(SDmanual)
}
