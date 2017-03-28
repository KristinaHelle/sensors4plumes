\name{optimiseSD_greedy}
\alias{optimiseSD_greedy}
\title{
Greedy optimisation algorithm
}
\description{
Runs greedy optimisation for a given Simulations object and cost function. The aim may be given as a number of sensors or as cost value, then the algorithm determines if sensors need to be added or deleted. Greedy optimisation adds (or respectively deletes) sensors one by one, always checking all possibilities and adding/deleting the sensor that yields minimal cost, given the already defined sensors; it stops when the aim is reached. In addition this algorithm can continue search when the greedy optimum is found by adding and deleting sensors in turns, this can help to get rid of redundant sensors and thus may improve the result.
The algorithm may take into account fix sensors and start from a given initial sensor set.
}
\usage{
optimiseSD_greedy(simulations, costFun,
  locationsAll = 1:nLocations(simulations), locationsFix = integer(0),
  locationsInitial = integer(0),
  aimCost = NA, aimNumber = NA,
  nameSave = NA, plot = FALSE, verbatim = FALSE,
  maxIterations = 100, swap = FALSE)
}

\arguments{
Some arguments are the same for all optimisation algorithms, they are marked by a \code{*}, for detail see \code{\link{optimiseSD}}
  \item{simulations}{\code{*}}
  \item{costFun}{\code{*}}
  \item{locationsAll}{\code{*}}
  \item{locationsFix}{\code{*}}
  \item{locationsInitial}{\code{*}}
  \item{aimCost}{\code{*}, ignored}
  \item{aimNumber}{\code{*}}
  \item{nameSave}{\code{*}}
  \item{plot}{not implemented yet}
  \item{verbatim}{\code{logical}, not implemented yet}
  \item{maxIterations}{maximal number of iterations, then it stops}
  \item{swap}{\code{logical} if algorithm continues, when greedy optimum is found, by adding and deleting sensors in turns}
}

\details{
In general the function is used within the wrapper \code{\link{optimiseSD}}. The parameters \code{maxIterations} and \code{swap} are specific to \code{optimiseSD_greedy}, they may be changed beforehand via \code{\link{replaceDefault}} with \code{type = "optimisationFun.optimiseSD"}; all other parameters are forwarded from \code{optimiseSD}.

If \code{aimCost} and \code{aimNumber} are given, \code{aimNumber} is ignored with a warning.
}
\value{
A list, the first two entries are common to all optimisation algorithms, they are marked with \code{*}, see \code{\link{optimiseSD}} for details.
\item{SD}{\code{*} best sampling designs}
\item{evaluation}{\code{*} cost and size of \code{SD}}
\item{report}{ a list of \cr
\code{SDs}a list of the sampling designs of all iterations \cr
\code{evalSDs} a \code{data.frame}: \code{cost} and \code{number} of sensors for the sampling designs \code{SDs}.}
}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
# optimisation function: is to be used inside of optimiseSD
# change parameters 'swap' and 'maxIterations'
optGr_20_true = replaceDefault(
  optimiseSD_greedy,
  newDefaults = list(
    maxIterations = 20,
    swap = TRUE)
)
# cost function
meanFun = function(x){mean(x, na.rm = TRUE)}
minDist = replaceDefault(
  spatialSpread, newDefaults = list(
    fun = minimalDistance,
    fun_R = meanFun
  ), type = "costFun.optimiseSD"
)[["fun"]]

# define possible, fix, and initial sensors
data(SimulationsSmall)
I = nLocations(SimulationsSmall)
set.seed(9345872)
locInit1 = sample.int(I, 2)
locKeep1 = sample(setdiff(1:I, locInit1), 2)
locAll1 = c(sample(setdiff(1:I, c(locInit1, locKeep1)), 4), locInit1)

# run optimisation
\dontrun{
## takes some time
optSD_greedy = optimiseSD(
  simulations = SimulationsSmall,
  costFun = minDist,
  optimisationFun = optGr_20_true[[1]],
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locInit1,
  aimNumber = 7
)

}
## this result is also in data(SDgreedy)
}
