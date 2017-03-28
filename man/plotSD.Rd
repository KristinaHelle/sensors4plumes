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
  locationsFix = integer(0), 
  locationsInitial = integer(0), 
  locationsAll = integer(0),
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
The function generates plots (main and point key missing) and returns a list of \code{trellis} objects (including main and keys).
}

\author{Kristina B. Helle, \email{kristina.helle@uni-muenster.de}}


\examples{
demo(radioactivePlumes_addProperties) 
data(SDssa)
SDs = list(SDssa[["SD"]][[1]][1,], SDssa$report[["SD_final"]])
singleDet = replaceDefault(singleDetection, 
  newDefaults = list(plot = TRUE), type = "costFun.optimiseSD")[[1]]
# separate maps
plotSD1a = plotSD(
  simulations = radioactivePlumes,
  SD = SDs,
  costMap = singleDet)
for (i in seq(along = SDs)){
  plot(plotSD1a[[i]])
}
  
# combined map (with customised parameters for fix and initial sensors)
I = nLocations(radioactivePlumes)
set.seed(22347287) # reconstruct initial and fix sensors of SDssa
locDel3 = sample.int(I, 5)
locKeep3 = sample(setdiff(1:I, locDel3), 10)
locAll3 = c(sample(setdiff(1:I, 
  c(locDel3, locKeep3)), 10), locDel3)

plotSD2 = plotSD(
  simulations = radioactivePlumes,
  SD = SDs,
  locationsAll = setdiff(1:nLocations(radioactivePlumes), c(locKeep3, locAll3)),
  locationsFix = locKeep3,
  locationsInitial = locDel3,
  costMap = singleDet,
  allIn1Plot = 2,
  col = c("green", "blue", "red", "white"),
  pch = c(1,2,4),
  col.regions = grey.colors,
  colorkey = FALSE)  


# combined map (with customised parameters for different SDs)   
plotSD3 = plotSD(
  simulations = radioactivePlumes,
  SD = SDs,
  locationsFix = locKeep3,
  locationsAll = setdiff(1:nLocations(radioactivePlumes), c(locKeep3, locAll3)), 
  allIn1Plot = 2,
  col = c(1, "white", "white", "white"),
  pch = c(1,20,0),
  col.regions = grey.colors,
  pointsKey = FALSE,
  pch.SDs = 2:5,
  cex.SDs = c(3.5,3,2.5,2),
  col.SDs = 4:8
  )  
}
