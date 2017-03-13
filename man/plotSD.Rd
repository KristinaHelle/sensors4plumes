\name{plotSD}
\alias{plotSD}

\title{
Plot sampling designs and related cost maps
}
\description{
Given a sampling design and a function to compute cost maps it generates the related map plot, based on \code{\link[sp]{spplot}}. If no cost map function is available, it may plot it with data from related \code{simulations}.

}
\usage{
plotSD(simulations, SD,
  locationsFix = integer(0), locationsInitial = integer(0), locationsAll = integer(0),
  costMap,
  zcol = 1, allIn1Plot = 0,
  pch = c(1, 20, 4), col = c("white", "white", "white", "white"), 
  pointsKey = TRUE, mainCost = TRUE,
  pch.SDs, col.SDs, cex.SDs, ...)
}

\arguments{
  \item{simulations}{ \code{Simulations} object}
  \item{SD}{ indices of \code{locations} of sensors, can be a matrix (rows refer to SDs) or a list}
  \item{locationsFix}{ indices of fix sensors, taken into account to compute cost map and plotted}
  \item{locationsInitial}{ indices of initial sensors, for plotting}
  \item{locationsAll}{ indices of all possible sensors, for plotting}
  \item{costMap}{ \code{function} to compute the values for the background map; needs the parameters \code{simulations} and \code{locations}; output must be a list with \code{"costLocations"} of length like the locations of \code{simulations}}
  \item{zcol}{ index or name: if no \code{costMap} given, this layer of the \code{locations} is used as background map}
  \item{allIn1Plot}{ if several SDs are given, they are by default (\code{allIn1Plot = 0}) plotted each on a map; if \code{allIn1Plot} is an integer, all SDs are combined in one plot with the cost map of this SD as background}
  \item{pch}{ point style of \code{SD}, \code{locationsFix}, \code{locationsInitial}; for \code{locationsAll} the \code{pch = "."} cannot be changed}
  \item{col}{ colors for \code{SD}, \code{locationsFix}, \code{locationsInitial}, \code{locationsAll}}
  \item{pointsKey}{logical, if key for points is to be printed}
  \item{mainCost}{logical, if cost of SD is to be used as main (this is only done if \code{costMap} returns a list item \code{"cost"})}
  \item{pch.SDs}{ point types for the different SDs, must have length that fits number of SDs (default: use \code{pch[1]} for all SDs but vary size)}
  \item{col.SDs}{ point colours for different SDs }
  \item{cex.SDs}{ point size for different SDs }
  \item{...}{ further parameters to be forwarded to \code{\link[sp]{spplot}}}
}

\value{
The function generates plots (main and point key missing) and returns a list of \code{\link[trellis]{trellis}} objects (including main and keys).
}

\author{Kristina B. Helle, \email{kristina.helle@uni-muenster.de}}


\examples{
data(radioactivePlumes_area)
# random SDs
SDs = list(
  sample.int(2500, 5),
  sample.int(2500, 10),
  sample.int(2500, 15),
  sample.int(2500, 20))

# prepare data
radioactivePlumes_area@values$detectable = calc(
  radioactivePlumes_area@values$maxdose,
  fun = function(x){x >= 1e-7})
radioactivePlumes_area@plumes$totalDose = 
  summaryPlumes(radioactivePlumes_area, fun = sum, values = "finaldose")[[2]]

# separate maps
plotSD1 = plotSD(
  simulations = radioactivePlumes_area,
  SD = SDs,
  costMap = singleDetection,
  mainCost = FALSE) # would not work, if all plots on one display 
c(plotSD1[[1]], plotSD1[[2]], plotSD1[[3]], plotSD1[[4]], layout = c(2,2))  # plot together

plotSD1a = plotSD(
  simulations = radioactivePlumes_area,
  SD = SDs,
  costMap = singleDetection)
for (i in seq(along = SDs)){
  plot(plotSD1a[[i]])
}
  
# combined map (with customised parameters for fix and initial sensors)
plotSD2 = plotSD(
  simulations = radioactivePlumes_area,
  SD = SDs,
  locationsFix = sample.int(2500, 30),
  locationsInitial = sample.int(2500, 40),
  locationsAll = sample.int(2500, 1000),
  costMap = singleDetection,
  allIn1Plot = 4,
  col = c("green", "blue", "red", "white"),
  pch = c(1,2,4),
  col.regions = grey.colors,
  colorkey = FALSE)  
plotSD2

# combined map (with customised parameters for different SDs)   
plotSD3 = plotSD(
  simulations = radioactivePlumes_area,
  SD = SDs,
  locationsFix = sample.int(2500, 30),
  locationsAll = sample.int(2500, 1000),,
  allIn1Plot = 4,
  col = c(1, "white", "white", "white"),
  pch = c(1,20,0),
  col.regions = grey.colors,
  pointsKey = FALSE,
  pch.SDs = 2:5,
  cex.SDs = c(3.5,3,2.5,2),
  col.SDs = 4:8
  )  



}