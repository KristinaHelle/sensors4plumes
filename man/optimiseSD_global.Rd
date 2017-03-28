\name{optimiseSD_global}
\alias{optimiseSD_global}
\title{
Derive one or all globally optimal sampling designs for plume detection
}
\description{
This optimisation works only for a special cost function: if the values are logical and indicating for each plume if it can be detected at any sensor location and cost is the fraction of plumes not detected by any sensor. In this case it can derive one or all global optima. It may take very long (hours to weeks) for cases of realistic size (hundreds of locations, more than 2 sensors). The algorithm is described in detail in the reference.
}
\usage{
optimiseSD_global(simulations, costFun = NA,
  locationsAll = 1:nLocations(simulations), locationsFix = integer(0),
  locationsInitial = integer(0),
  aimCost = NA, aimNumber = NA,
  nameSave = NA, plot = FALSE, verbatim = FALSE,
  detectable = 1, maxIterations = NA,
  findAllOptima = FALSE, findSensorNumber = FALSE
  )
}

\arguments{
Some arguments are the same for all optimisation algorithms, they are marked by a \code{*}, for detail see \code{\link{optimiseSD}}
  \item{simulations}{\code{*}}
  \item{costFun}{\code{*} ignored, as optimisation can only optimise plume detection}
  \item{locationsAll}{\code{*}}
  \item{locationsFix}{\code{*}}
  \item{locationsInitial}{\code{*} ignored, only if \code{aimNumber} missing it is defined by length of \code{locationsInitial}}
  \item{aimCost}{\code{*}, ignored}
  \item{aimNumber}{\code{*}}
  \item{nameSave}{\code{*}}
  \item{plot}{not implemented yet}
  \item{verbatim}{if more output is to be written to the console}
  \item{detectable}{which layer of the values of \code{simulations} to use, character or integer}
  \item{maxIterations}{stop when this number of sensor sets was tested (if \code{findAllOptima = TRUE} this is applied to each search i.e. to the one to find the maximal number of detectable plumes and to the one to find all optima); if \code{maxIterations} is below the number of sensor sets to actually be tested, the result may not be the global optimum.}
  \item{findAllOptima}{logical, if all optima to be found - requires to internally run search twice.}
  \item{findSensorNumber}{logical, only relevant if \code{findAllOptima = FALSE} and all plumes are detected; in this case \code{findSensorNumber = TRUE} calls further searches to determine the minimal number of sensors required to detect all plumes. If \code{findAllOptima = TRUE} these searches are run in any case.}
}

\details{
In general the function is used within the wrapper \code{\link{optimiseSD}}. The non-\code{*}-parameters are specific to \code{optimiseSD_global}, they may be changed beforehand via \code{\link{replaceDefault}} with \code{type = "optimisationFun.optimiseSD"}; all other parameters are forwarded from \code{optimiseSD}.


If \code{aimCost} is given, the algorithm stops when the aim is reached by the best sampling design. In this case, nothing is returned, but the populations are saved (at \code{nameSave} and if nothing indicated at \code{"opt_genetic.Rdata"}), see examples how to extract the sampling designs in this case.
}

\value{
A list, the first two entries are common to all optimisation algorithms, they are marked with \code{*}, see \code{\link{optimiseSD}} for details.
\item{SD}{\code{*} best sampling designs}
\item{evaluation}{\code{*} cost and size of \code{SD}}
\item{report}{a list with more details about the optimisation\cr
\code{detectable} the matrix of 0 and 1 used in the optimisation after selecting only locations in \code{locationsAll}, deleting locations in \code{locationsFix} and containing only plumes that cannot be detected at \code{locationsFix} but at \code{locationsAll}\cr
\code{first} results of the first run to find the detection limit, \code{lowerLimit} is this limit plus 1\cr
\code{all} results of the run to find all optimal sampling designs
}
}
\references{
K.B. Helle, E. Pebesma (2015). Optimising sampling designs for the maximum coverage problem of plume detection. Spatial Statistics (13), 21-44.
}

\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
# prepare data
# test case from reference
detectionMatrix = matrix(c(0,1,0,0,0,1,
                          1,0,1,0,0,0,
                          0,1,1,1,0,0,
                          1,0,0,0,0,0,
                          0,0,0,1,1,0,
                          1,1,0,0,0,0,
                          0,0,1,0,0,1),
                        nrow = 7, ncol = 6, byrow = TRUE)
data(SimulationsSmall)
completeExample = Simulations(
  locations = SimulationsSmall@locations[1:7,],
  plumes = SimulationsSmall@plumes[1:6,],
  values = raster(x = detectionMatrix,
             xmn = -90, xmx = 90, ymn = -90, ymx = 90,
             crs = "+init=epsg:4326"))

# the function is to be used inside of optimiseSD
# change algorithm specific parameters 'detectable', 'findAllOptima'
optSDglobal = replaceDefault(
  fun = optimiseSD_global,
  newDefaults = list(
    detectable = 1,
    findAllOptima = TRUE
    ),
  type = "optimisationFun.optimiseSD"
)

# run optimisation
\dontrun{
# takes some (little) time
optSD_global = optimiseSD(
  simulations = completeExample,
  aimNumber = 3,
  costFun = NA,
  optimisationFun = optSDglobal[[1]],
  nameSave = NA
)
}
## this result is also in data(SDglobal)

data(SDglobal)
# visualise result
## which plumes are detected?
detAtSensors = matrix(completeExample@values[SDglobal$SD[[1]][1,],],
  byrow = TRUE, nrow = 3)
undetedtedPlumes = ! apply(FUN = any, X = detAtSensors, MARGIN = 2)
## cost is fraction of undetected plumes
sum(undetedtedPlumes)/nPlumes(completeExample) == SDglobal$cost

## map of undetected plumes
detectablePlumesLoc = matrix(getValues(completeExample@values),
  byrow = TRUE, ncol = nPlumes(completeExample))
completeExample@locations@data$detectable = apply(FUN = sum, X = detectablePlumesLoc,
	MARGIN = 1)

SDopt = as(subsetSDF(completeExample@locations,
                   locations =  SDglobal$SD[[1]][1,]),
         "SpatialPointsDataFrame")

spplot(completeExample@locations, zcol = "detectable",
  sp.layout = list(list("sp.points",
  SpatialPoints(coords = coordinates(SDopt), proj4string = CRS(proj4string(SDopt))),
                        col = 3, cex = 2, lwd = 1.5)))
}
