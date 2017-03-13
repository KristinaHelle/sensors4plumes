\name{optimiseSD_genetic}
\alias{optimiseSD_genetic}
\title{
Optimisation of sensor locations by a binary genetic algorithm 
}
\description{
Optimises sensors for a given \code{Simulations} object and \code{costFunction}. It starts from a population of sensor sets that can be provided and else are random. By testing them and recombining good locations, a new improved population is generated. This process continues until the aim or a stopping criterion is reached. It applies the binary genetic algorithm \code{\link[genalg]{rbga.bin}}.   
}
\usage{
optimiseSD_genetic(simulations, costFun, 
  locationsAll = 1:nLocations(simulations), locationsFix = integer(0), locationsInitial = integer(0),
  aimCost = NA, aimNumber = NULL,
  nameSave = NA, verbatim = FALSE,
  evalFunc = numberPenalty, 
  popSize = 200, iters = 100, 
  mutationChance = NA, elitism = NA
)
}

\arguments{
Some arguments are the same for all optimisation algorithms, they are marked by a \code{*}, for detail see \code{\link{optimiseSD}}. Others are directly forwarded to \code{\link[genalg]{rbga.bin}}, they are marked by \code{**}
  \item{simulations}{\code{*}}
  \item{costFun}{\code{*}}
  \item{locationsAll}{\code{*}}
  \item{locationsFix}{\code{*}}
  \item{locationsInitial}{\code{*}; here it may also be a matrix, each row indicates a set of sensors}
  \item{aimCost}{\code{*}}
  \item{aimNumber}{\code{*}; is ignored if \code{locationsInitial} given, else it defines the number of sensors in the sets of the initial population}
  \item{nameSave}{\code{*}}
  \item{verbatim}{if more output is to be written to the console}
  \item{evalFunc}{\code{**}; this function has to take a chromosome and return its cost (including fix sensors); it can be modified by \code{\link{replaceDefault}} with \code{type = "evalFunc.rbga.bin"}}
  \item{popSize}{\code{**}} 
  \item{iters}{\code{**}} 
  \item{mutationChance}{\code{**}}
  \item{elitism}{\code{**}}
}

\details{
In general the function is used within the wrapper \code{\link{optimiseSD}}. The \code{**}-parameters are specific to \code{optimiseSD_greedy}, they may be changed beforehand via \code{\link{replaceDefault}} with \code{type = "optimisationFun.optimiseSD"}; all other parameters are forwarded from \code{optimiseSD}. 
If \code{aimCost} is given, the algorithm stops when the aim is reached by the best sampling design. In this case, nothing is returned, but the populations are saved (at \code{saveName} and if nothing indicated at \code{"opt_genetic.Rdata"}), see examples how to extract the sampling designs in this case.
}

\value{
A list
\item{SD}{a list of all best designs}
\item{cost}{cost of all \code{SD}}
\item{report}{the generated \code{"rbga"} object as returned by \code{\link{rbga.bin}}}
}
\author{ kristina.helle@uni-muenster.de  }

\examples{
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
OptSD_gen1 = optimiseSD(
  simulations = radioactivePlumes_local, 
  costFun = singleDetection,
  optimisationFun = optimiseSD_genetic_1,
  locationsAll = locAll_l, 
  locationsFix = locKeep_l, 
  locationsInitial = locInit_l[1,],
  nameSave = NA, # needed, as else NULL is forwarded and parameter destroyed 
  aimCost = NA # needed, as else NULL is forwarded and parameter destroyed 
)

## case with 'aimCost'
optimiseSD_genetic_2 = replaceDefault(
    optimiseSD_genetic, newDefaults = list(
      popSize = 20, 
      iters = 20
    ),
    type = "optimisationFun.optimiseSD")[[1]]

set.seed(07021916)
OptSD_gen2 = optimiseSD(
  simulations = radioactivePlumes_local, 
  costFun = singleDetection,
  optimisationFun = optimiseSD_genetic_2,
  locationsAll = locAll_l, 
  locationsFix = locKeep_l, 
  locationsInitial = locInit_l[1,],
  aimCost = 0.05,
  nameSave = "optSD_genetic_2"
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
    costs[j,i] = evalFun(populations[i,,j])
    SDs[[j]][[i]] = c(locAll_l[populations[i,,j] == 1], locKeep_l)
    if (sum(populations[i,,j]) > 5){
      costs[j,i] = 2  
    }
    print(paste(j,i))
  }  
}
apply(FUN = min, X = costs, MARGIN = 1) # best cost in each iteration
SDs[[finalIteration]][[which(costs[finalIteration,] == min(costs[finalIteration,]))]] # best sampling design
}
