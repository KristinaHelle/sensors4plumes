\name{optimiseSD}
\alias{optimiseSD}
\title{
(Spatial) optimisation of sampling designs
}
\description{
This function optimises a sampling design to achieve minimal cost for a given cost function by applying a given optimisation function. It may take into account fix locations or constrain search to a subset of locations.
}
\usage{
optimiseSD(simulations, costFun,
  locationsAll = 1:nLocations(simulations), locationsFix = integer(0), locationsInitial = integer(0),
  aimCost, aimNumber, 
  optimisationFun,
  nameSave ="optimiseSD", plot = FALSE, verbatim = FALSE, ...)
}

\arguments{
  \item{simulations}{ \code{\link{Simulations} object}}
  \item{costFun}{cost function, must have parameters \code{simulations} and \code{locations} (may be prepared by \code{replaceDefault} with \code{type = "costFun.optimiseSD"}); it must return cost as a single value or a list where this is the first entry}
  \item{locationsAll}{indices of the \code{locations} that are considered possible sensor locations, by default these are all locations}
  \item{locationsFix}{indices of \code{locations} with fix sensors -- to be considered when computing cost, by default this is empty}
  \item{locationsInitial}{indices of \code{locations} where sensors are initially; by default this is empty; if not indicated \code{aimNumber} random locations may be used, depending on the }
  \item{aimCost}{limit value of the \code{costFun}: the result has to fall below}
  \item{aimNumber}{if no \code{locationsInitial} given, it starts from this number of random locations}
  \item{optimisationFun}{\code{function} that executes an optimisation algorithm; must have parameters \code{simulations, costFun, locationsAll, locationsFix, locationsInitial, aimCost, aimNumber, nameSave} (use \code{type = "optimisationFun.optimiseSD"} in \code{\link{replaceDefault}}). Output must be a list with \code{SD} (vector, matrix, or list of one or several SDs) and \code{cost} (belonging to the SD); it may return other values in \code{report}}
  \item{nameSave}{\code{character} path and name (without suffix!) where to save intermediate results (sampling design of each iteration)}
  \item{plot}{if iteration is plotted (currently disabled)}
  \item{verbatim}{print intermediate results and keep sampling desings and cost of all iterations}
  \item{...}{further parameters, currently unused}
}

%\details{
%}
%\value{
%A list
%\item{SD}{final sampling desing}
%\item{SD_best}{best sampling design ever found}
%\item{cost}{cost of \code{SD}}
%\item{cost_best}{cost of \code{SD_best}}
%\item{report}{\code{data.frame} with columns \code{cost} and \code{accepted}, rows refer to iterations.}

%If \code{verbatim = TRUE} more information is returned:
%The \code{report} has additional columns: \code{costTest} (cost of the design tested in this iteration, no matter if it was accepted; \code{cost} refers to the currently accepted design); \code{costBest} (best cost until then); \code{chi} (random value that was used for the acceptance decision).

%In addition it returns a matrix each for \code{SDs}, \code{SDs_test}, \code{SDs_test} with the current, tested and best known sampling designs of all iterations.

%If \code{acceptanceMethod = "vanGroenigen"} it returns the cooling parameters \code{cooling} and \code{startAcceptance} as they may be computed inside the algorithm.

%It gives the iterations with \code{jumpBack} i.e. when the current sampling design was replaced by the known best.
%}
\author{ kristina.helle@uni-muenster.de  }

\examples{
# the function is to be used inside of optimiseSD
# change parameters 
optimSD_ssa1 = replaceDefault(
  optimiseSD_ssa,  newDefaults = list(
    start_acc_vG = 0.1,
    aimCost = 0, 
    verbatim = TRUE,
    maxIterations = 3000,
    maxStableIterations = 500,
    maxIterationsJumpBack = 200   
  ),
  type = "optimisationFun.optimiseSD")[[1]] 

# define possible, fix, and initial sensors
data(radioactivePlumes_area)
locDel3 = sample.int(nLocations(radioactivePlumes_area), 5)
locKeep3 = sample(setdiff(1:nLocations(radioactivePlumes_area), locDel3), 100)
locAll3 = c(sample(setdiff(1:nLocations(radioactivePlumes_area), 
  c(locDel3, locKeep3)), 10), locDel3)

# prepare data
threshold = 1e-7
radioactivePlumes_area@values$detectable = calc(
  radioactivePlumes_area@values$maxdose,
  fun = function(x){x >= threshold})
radioactivePlumes_area@plumes$nDetectable = 
  summaryPlumes(radioactivePlumes_area, fun = sum, values = "detectable")[[2]]  
costInitial1 = multipleDetection(simulations = radioactivePlumes_area, 
  locations = c(locKeep3, locDel3))

# run optimisation
optSSA1 = optimiseSD(
  simulations = radioactivePlumes_area,
  costFun = multipleDetection,
  locationsAll = setdiff(1:nLocations(radioactivePlumes_area), c(locKeep3, locAll3)),                       
  locationsFix = locKeep3,
  locationsInitial = locDel3,
  aimCost = 0.05 * costInitial1[[1]], 
  optimisationFun = optimSD_ssa1
)  

# visualise 
## cost curve
plot(optSSA1$report$costTest, pch = ".")
lines(optSSA1$report$costBest, col = 4)
lines(optSSA1$report$cost, col = 3)
points(optSSA1$report$costTest, pch = ".", col = optSSA1$report$accepted + 2)
## designs
allSDs = subsetSDF(radioactivePlumes_area@locations,
  locations = c(locKeep3, locDel3, optSSA1$SD[1:5], unique(as.integer(optSSA1$SDs_test))))
AllSDs = as(allSDs, "SpatialPointsDataFrame")  
radioactivePlumes_area@locations@data$nDetectable = 
  summaryLocations(radioactivePlumes_area, fun = sum, values = "detectable")[[2]] 
spplot(radioactivePlumes_area@locations, zcol = "nDetectable",
  sp.layout = list("sp.points", AllSDs, 
  pch = c(rep(4, length(locKeep3)), 
          rep(1, length(locDel3)), 
          rep(20, 5), rep(46, 1450)), 
  col = "green"))

}
