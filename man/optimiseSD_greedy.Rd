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
  locationsAll, locationsFix, locationsInitial, 
  aimCost, aimNumber,
  maxIterations = 100, swap = FALSE, 
  nameSave)
}

\arguments{
  \item{simulations}{ \code{\link{Simulations} object}}
  \item{costFun}{cost function, must have parameters \code{simulations} and \code{locations} (may be prepared by \code{replaceDefault} with \code{type = "costFun.optimiseSD"}); it must return cost as a single value or a list where this is the first entry}
  \item{locationsAll}{indices of the \code{locations} that are considered possible sensor locations}
  \item{locationsFix}{indices of \code{locations} with fix sensors -- to be considered when computing cost}
  \item{locationsInitial}{indices of \code{locations} where sensors are initially}
  \item{aimCost}{limit value of the \code{costFun}: the result has to fall below}
  \item{aimNumber}{desired number of sensors; the result has to have no more sensors (\code{locationsFix} do not count)}
  \item{maxIterations}{maximal number of iterations, then it stops}
  \item{swap}{\code{logical} if algorithm continues, when greedy optimum is found, by adding and deleting sensors in turns}
  \item{nameSave}{\code{character} path and name (without suffix!) where to save intermediate results (sampling design of each iteration)}
}

\details{
In general the function is used within the wrapper \code{\link{optimiseSD}}. The parameters \code{maxIterations} and \code{swap} are specific to \code{optimiseSD_greedy}, they may be changed beforehand via \code{\link{replaceDefault}} with \code{type = "optimisationFun.optimiseSD"}; all other parameters are forwarded from \code{optimiseSD}. 

If \code{aimCost} and \code{aimNumber} are given, \code{aimNumber} is ignored with a warning.
}
\value{
A list
\item{SDs}{a list of the sampling designs of all iterations}
\item{finalSDwhich}{one or more integers indicating the iteration of the best sampling design(s): they fulfill the given aim and have minimal number of sensors/minimal cost. If the aim cannot be reached (\code{maxIterations} reached; all possible sensors added and still too few or too high cost; all but fix sensors deleted but still too many / too low cost) it may be empty; if cost too low or too few sensors the sampling design with highest cost/most sensors is returned.}
\item{evalSDs}{\code{data.frame}: \code{cost} and \code{number} of sensors for the sampling design of each iteration.}
}
\author{ kristina.helle@uni-muenster.de  }

\examples{
# the function is to be used inside of optimiseSD
# change parameters 'swap' and 'maxIterations'
optGr_20_true = replaceDefault(
  optimiseSD_greedy,
  newDefaults = list(
    maxIterations = 20,
    swap = TRUE)
)
# define possible, fix, and initial sensors
data(SimulationsSmall)
locInit1 = sample.int(nLocations(SimulationsSmall), 2)
locKeep1 = sample(setdiff(1:nLocations(SimulationsSmall), locInit1), 2)
locAll1 = c(sample(setdiff(1:nLocations(SimulationsSmall), c(locInit1, locKeep1)), 4), locInit1)

# run optimisation
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
