\name{optimiseSD_genetic}
\alias{optimiseSD_genetic}
\alias{numberPenalty}

\title{
Optimisation of sensor locations by a binary genetic algorithm
}
\description{
\code{optimiseSD_genetic} optimises sensors for a given \code{Simulations} object and \code{costFunction}. It starts from a population of sensor sets that can be provided and else are random. By testing them and recombining good locations, a new improved population is generated. This process continues until the aim or a stopping criterion is reached. It applies the binary genetic algorithm \code{\link[genalg]{rbga.bin}}.\cr
\code{numberPenalty} is the default evaluation function, it returns the penalty if \code{aimNumber} is exceeded, else the cost computed by \code{costFun}.
}
\usage{
optimiseSD_genetic(simulations, costFun,
  locationsAll = 1:nLocations(simulations), locationsFix = integer(0),
  locationsInitial = integer(0),
  aimCost = NA, aimNumber = NA,
  nameSave = NA, plot = FALSE, verbatim = FALSE,
  evalFunc = numberPenalty, 
  popSize = 200, iters = 100, 
  mutationChance = NA, elitism = NA)
numberPenalty(chromosome, simulations, costFun, 
  locationsAll, locationsFix, 
  aimNumberNF, penalty = 2)
}

\arguments{
Some arguments are the same for all optimisation algorithms, they are marked by a \code{*}, for detail see \code{\link{optimiseSD}}. Others are directly forwarded to \code{\link[genalg]{rbga.bin}}, they are marked by \code{**}
  \item{simulations}{\code{*}}
  \item{costFun}{\code{*}}
  \item{locationsAll}{\code{*}}
  \item{locationsFix}{\code{*}}
  \item{locationsInitial}{\code{*}; here it may also be a matrix, each row indicates a set of sensors}
  \item{aimCost}{\code{*}}
  \item{aimNumber}{\code{*}; if \code{locationsInitial} given it must be \code{NULL} and is ignored, else it defines the number of sensors in the sets of the initial population}
    \item{aimNumberNF}{\code{aimNumber} minus length of \code{locationsFix} - to be used in evaluation}
  \item{nameSave}{\code{*}}
  \item{plot}{not yet implemented}
  \item{verbatim}{if more output is to be written to the console}
  \item{evalFunc}{\code{**}; this function has to take a chromosome and return its cost (including fix sensors); it can be modified by \code{\link{replaceDefault}} with \code{type = "evalFunc.rbga.bin"}}
  \item{popSize}{\code{**}}
  \item{iters}{\code{**}}
  \item{mutationChance}{\code{**}}
  \item{elitism}{\code{**}}
  \item{chromosome}{a vector of 0 and 1 of the same length as \code{locationsAll} indicating sensors at the positions of 1.}
  \item{penalty}{\code{numeric} to be returned if number of 1s in \code{chromosome} exceeds \code{aimNumber}.}
}

\details{
In general the function is used within the wrapper \code{\link{optimiseSD}}. The \code{**}-parameters are specific to \code{optimiseSD_greedy}, they may be changed beforehand via \code{\link{replaceDefault}} with \code{type = "optimisationFun.optimiseSD"}; all other parameters are forwarded from \code{optimiseSD}.
If \code{aimCost} is given, the algorithm stops when the aim is reached by the best sampling design. In this case, nothing is returned, but the populations are saved (at \code{nameSave} and if nothing indicated at \code{"opt_genetic.Rdata"}), see examples how to extract the sampling designs in this case.
}

\value{
A list, the first two entries are common to all optimisation algorithms, they are marked with \code{*}, see \code{\link{optimiseSD}} for details.
\item{SD}{\code{*} best sampling designs}
\item{evaluation}{\code{*} cost and size of \code{SD}}
\item{report}{the generated \code{"rbga"} object as returned by \code{\link{rbga.bin}}}
}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
# load data and compute required properties
demo(radioactivePlumes_addProperties)

# initial, fix, possible sensor locations
I = nLocations(radioactivePlumes)
set.seed(1000)
locInit_l = t(replicate(10, sample.int(I, 5))) 
locKeep_l = sample(setdiff(1:I, locInit_l), 2)
locAll_l = c(sample(setdiff(1:I, c(locInit_l, locKeep_l)),
                          round(I - (5 + 2)) * 0.5), locInit_l)

# the function is to be used inside of optimiseSD
# change algorithm specific parameters 'evalFunc', 'popSize', 'iters'
optimiseSD_genetic_1 = replaceDefault(
  optimiseSD_genetic, newDefaults = list(
    evalFunc = numberPenalty,
    popSize = 20,
    iters = 10
  ),
  type = "optimisationFun.optimiseSD")[[1]]

# run optimisation
\dontrun{
## takes some time
OptSD_gen1 = optimiseSD(
  simulations = radioactivePlumes,
  costFun = singleDetection,
  optimisationFun = optimiseSD_genetic_1,
  locationsAll = locAll_l,
  locationsFix = locKeep_l,
  locationsInitial = locInit_l[1,])
}
## this result is also in data(SDgenetic)

\dontrun{
# result is not returned but saved to file
## case with 'aimCost'
optimiseSD_genetic_2 = replaceDefault(
    optimiseSD_genetic, newDefaults = list(
      evalFunc = numberPenalty,
      popSize = 20,
      iters = 20
    ),
    type = "optimisationFun.optimiseSD")[[1]]

set.seed(07021916)
OptSD_gen2 = optimiseSD(
  simulations = radioactivePlumes,
  costFun = singleDetection,
  optimisationFun = optimiseSD_genetic_2,
  locationsAll = locAll_l,
  locationsFix = locKeep_l,
  locationsInitial = locInit_l[1,],
  aimCost = 0.05,
  nameSave = "optSD_genetic_2" # result is saved to file here
)
# OptSD_gen2 not found
## load generated populations and extract sampling designs and cost
load("optSD_genetic_2.Rdata")
finalIteration = sum(!is.na(populations[1,1,]))
costs = matrix(0, nrow = finalIteration, ncol = 20)
SDs = list()
for (j in 1:finalIteration){
  SDs[[j]] = list()
  for (i in 1:20){
    costs[j,i] = numberPenalty(populations[i,,j],
      simulations = radioactivePlumes,
      costFun = singleDetection,
      locationsAll = locAll_l,
      locationsFix = locKeep_l,
      aimNumber = length(locInit_l[1,]))
    SDs[[j]][[i]] = c(locAll_l[populations[i,,j] == 1], locKeep_l)
    if (sum(populations[i,,j]) > 5){
      costs[j,i] = 2
    }
    print(paste(j,i))
  }
}
apply(FUN = min, X = costs, MARGIN = 1) # best cost in each iteration
# best sampling design:
SDs[[finalIteration]][[which(costs[finalIteration,] == min(costs[finalIteration,]))]]
}
}
