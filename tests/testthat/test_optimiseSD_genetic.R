################################################################
# test optimiseSD_genetic                                      #
################################################################
# evalFun: numberPenalty
# rgba2SD 
data(radioactivePlumes_local)
locInit_l = t(replicate(10, sample.int(nLocations(radioactivePlumes_local), 5))) # samples in rows (replicate: samples in columns)
locKeep_l = sample(setdiff(1:nLocations(radioactivePlumes_local), locInit_l), 10)
locAll_l = c(sample(setdiff(1:nLocations(radioactivePlumes_local), c(locInit_l, locKeep_l)), 
                          nLocations(radioactivePlumes_local) - 60 - 1000), locInit_l)
threshold = 1e-7
radioactivePlumes_local@values$detectable = calc(
     radioactivePlumes_local@values$maxdose,
     fun = function(x){x >= threshold})
radioactivePlumes_local@plumes$totalDose = 
     summaryPlumes(radioactivePlumes_local, fun = sum, values = "finaldose")[[2]]
radioactivePlumes_local@plumes$nDetectable = 
     summaryPlumes(radioactivePlumes_local, fun = sum, values = "detectable")[[2]]


# ------------------- numberPenalty ------------------------ #
cost0_singleDet_l = singleDetection(simulations = radioactivePlumes_local, locations = integer(0))
chromosome1 = rbinom(length(locAll_l), 1, 5/length(locAll_l))
costNP1 = numberPenalty(
  simulations = radioactivePlumes_local,
  costFun = singleDetection,
  chromosome = chromosome1,
  locationsAll = locAll_l,
  locationsFix = locKeep_l,
  aimNumber = 15,
  penalty = 2 * cost0_singleDet_l[[1]],
  plot = FALSE)

if (sum(chromosome1) > 5){
  expect_equal(
    costNP1[[1]],
    2 * cost0_singleDet_l[[1]]
  )  
} else {
  expect_equal(
    costNP1[[1]],
    singleDetection(
      simulations = radioactivePlumes_local,
      locations = c(locKeep_l,  locAll_l[chromosome1 == 1])
      )[[1]]
  )  
}

#------------ optimiseSD_genetic ------------------------- #
# set seed and determine if result is as expected (do not evaluate what happens inside rbga.bin)
set.seed(09112006)
optSD_gen1 = optimiseSD_genetic(
  simulations = radioactivePlumes_local, 
  costFun = singleDetection, 
  evalFun = numberPenalty,
  locationsAll = locAll_l, 
  locationsFix = locKeep_l, 
  locationsInitial = locInit_l[1,], # vector
  popSize = 20, 
  iters = 10,
  verbatim = TRUE
)

## all input forwarded as expected
evaluationFunction_singleDet = replaceDefault(
  fun = numberPenalty,
  newDefaults = list(
    simulations = radioactivePlumes_local,
    costFun = singleDetection,
    locationsAll = locAll_l,
    locationsFix = locKeep_l,
    aimNumber = ncol(locInit_l)
  )
)[[1]]
suggestions_locInit_l = matrix(0, nrow = nrow(locInit_l), ncol = length(locAll_l))
for (i in 1:nrow(locInit_l)){
  suggestions_locInit_l[i,is.element(locAll_l, locInit_l[i,])] = 1   
}

set.seed(09112006)
optSD_Gen1 = rbga.bin(
  size = length(locAll_l), 
  evalFun = evaluationFunction_singleDet,
  popSize = 20, 
  iters = 10,
  suggestions = suggestions_locInit_l[1,,drop = FALSE],
  zeroToOneRatio = round(length(locAll_l)/ncol(locInit_l)),
)

expect_equal(
  optSD_gen1$report,
  optSD_Gen1
)

set.seed(12022014)
optSD_gen2 = optimiseSD_genetic(
  simulations = radioactivePlumes_local, 
  costFun = multipleDetection, 
  evalFun = numberPenalty,
  locationsAll = locAll_l, 
  locationsFix = locKeep_l, 
  locationsInitial = locInit_l, # matrix
  popSize = 20, 
  iters = 10,
  mutationChance = 0.01,
  elitism = 1
)

# size fits locationsAll 
expect_equal(
  optSD_gen2$report$size,
  length(locAll_l) # may be wrong if locAll overlaps with locFix or does not cover all locInitial
)
# cost
i = sample.int(optSD_gen2$report$popSize, 1)
## computed by evalFun (including costFun); for correct use of locations see test of numberPenalty
## aimNumber is length of locationsInitial 
expect_equal(
  optSD_gen2$report$evaluations[i],
  numberPenalty(
    simulations = radioactivePlumes_local,
    costFun = multipleDetection,
    chromosome = optSD_gen2$report$population[i,],
    locationsAll = locAll_l,
    locationsFix = locKeep_l,
    aimNumber = length(locInit_l),
    penalty = 2,
    plot = FALSE)
)

## all input forwarded as expected
evaluationFunction_multipleDet = replaceDefault(
  fun = numberPenalty,
  newDefaults = list(
    simulations = radioactivePlumes_local,
    costFun = multipleDetection,
    locationsAll = locAll_l,
    locationsFix = locKeep_l,
    aimNumber = ncol(locInit_l)
  )
)[[1]]

set.seed(12022014)
optSD_Gen2 = rbga.bin(
  size = length(locAll_l), 
  evalFun = evaluationFunction_multipleDet,
  popSize = 20, 
  iters = 10,
  suggestions = suggestions_locInit_l,
  zeroToOneRatio = round(length(locAll_l)/ncol(locInit_l)),
  mutationChance = 0.01,
  elitism = 1
)

expect_equal(
  optSD_gen2$report,
  optSD_Gen2
)

set.seed(02082007)
optSD_gen3 = optimiseSD_genetic(
  simulations = radioactivePlumes_local, 
  costFun = singleDetection, 
  evalFun = numberPenalty,
  locationsAll = locAll_l, 
  locationsFix = locKeep_l, 
  # no locationsInitial
  aimNumber = 5,
  popSize = 20, 
  iters = 10,
  verbatim = TRUE
)

## all input forwarded as expected
set.seed(02082007)
locInit = t(replicate(19, sample(locAll_l, 5)))
newLocations = cleanIndices(
  locationsTotal = 1:nLocations(radioactivePlumes_local),
  locationsAll = locAll_l,
  locationsFix = locKeep_l,
  locationsInitial = locInit)
suggestionsInit = matrix(0, nrow = nrow(newLocations$locationsInitial), 
                         ncol = length(newLocations$locationsAll))
for (i in 1:dim(suggestionsInit)[1]){
  suggestionsInit[i, is.element(newLocations$locationsAll, newLocations$locationsInitial[i,])] = 1
}
optSD_Gen3 = rbga.bin(
  size = length(newLocations$locationsAll), 
  evalFun = evaluationFunction_singleDet,
  popSize = 20, 
  iters = 10,
  suggestions = suggestionsInit,
  zeroToOneRatio = round(length(newLocations$locationsAll)/ncol(newLocations$locationsInitial)),
)

expect_equal(
  optSD_gen3$report,
  optSD_Gen3
)

# -------------------- monitor ----------------------------- #
## save if nameSave given
file.remove("gen1.Rdata")
optSD_gen_m1 = optimiseSD_genetic(
  simulations = radioactivePlumes_local, 
  costFun = singleDetection, 
  evalFun = numberPenalty,
  locationsAll = locAll_l, 
  locationsFix = locKeep_l, 
  locationsInitial = newLocations$locationsInitial,
  nameSave = "gen1",
  popSize = 20, 
  iters = 10,
  verbatim = TRUE
)
load("gen1.Rdata")
expect_equal(# first entry of populations contains suggestions
    populations[1:19,,1],
    suggestionsInit    
)
expect_equal(# last entry of populations equals final population
  populations[,,10],
  optSD_gen_m1$population
)

## stop when aimCost reached (save automatically)
set.seed(07021916)
optSD_gen_m2 = optimiseSD_genetic(
  simulations = radioactivePlumes_local, 
  costFun = singleDetection, 
  evalFun = numberPenalty,
  locationsAll = locAll_l, 
  locationsFix = locKeep_l, 
  locationsInitial = locInit_l[1,],
  aimCost = 0.05, 
  popSize = 20, 
  iters = 20
)
load("opt_genetic.Rdata")
# extract SD and compute cost
evalFun = replaceDefault(
  fun = numberPenalty,
  newDefaults = list(
    simulations = radioactivePlumes_local,
    costFun = singleDetection,
    locationsAll = locAll_l,
    locationsFix = locKeep_l,
    aimNumber = 5
    #plot = plot
  ), type = "evalFunc.rbga.bin"
)[[1]]

finalIteration = sum(!is.na(populations[1,1,]))
costs = matrix(0, nrow = finalIteration, ncol = 20)
SDs = list()
for (j in 1:finalIteration){
  SDs[[j]] = list()
  for (i in 1:20){
    costs[j,i] = evalFun(populations[i,,j])
    if (sum(populations[i,,j]) > 5){
      costs[j,i] = 2  
    }
    print(paste(j,i))
  }  
}
apply(FUN = min, X = costs, MARGIN = 1)

expect_less_than(# best cost of last iteration is below limit
  min(costs[finalIteration,]),
  0.05
)
expect_more_than(# before last iteration, cost never below limit
  min(costs[1:(finalIteration-1),]),
  0.05
)

## (plot)
## (user-defined monitoring)

# ------------- rbga2SD ----------------------------------- #

# transforming back suggestions returns initial locations
# cost of any population member is as expected

# SDs are the best SDs in the population
costs_gen1 = numeric(20)
SDs_gen1 = list()
for (i in 1:20){
  SDs_gen1[[i]] = c(locAll_l[optSD_gen1$report$population[i,] == 1], locKeep_l)
  costs_gen1[i] = singleDetection(simulations = radioactivePlumes_local,
                             locations = SDs_gen1[[i]])[[1]]
  if (sum(optSD_gen1$report$population[i,]) > 5){
    costs_gen1[i] = 2  
  }
}  
best_gen1 = which(costs_gen1 == min(costs_gen1))
uniqueBest_gen1 = unique(optSD_gen1$report$population[best_gen1,])
expect_equal(# evaluation is cost of SDs
  optSD_gen1$report$evaluation,
  costs_gen1
)
expect_equal(# same number of best designs
  dim(uniqueBest_gen1)[1],
  length(optSD_gen1$SD)
)
for (i in 1:dim(uniqueBest_gen1)[1]){
  expect_equal(# same best designs
    sort(c(locAll_l[uniqueBest_gen1[i,] == 1], locKeep_l)),
    sort(optSD_gen1$SD[[i]])
  )  
}


costs_gen2 = numeric(20)
SDs_gen2 = list()
for (i in 1:20){
  SDs_gen2[[i]] = c(locAll_l[optSD_gen2$report$population[i,] == 1], locKeep_l)
  costs_gen2[i] = multipleDetection(simulations = radioactivePlumes_local,
                                  locations = SDs_gen2[[i]])[[1]]
  if (sum(optSD_gen2$report$population[i,]) > 5){
    costs_gen2[i] = 2  
  }
}
best_gen2 = which(costs_gen2 == min(costs_gen2))
uniqueBest_gen2 = unique(optSD_gen2$report$population[best_gen2,])

expect_equal(# evaluation is cost of SDs
  optSD_gen2$report$evaluation,
  costs_gen2
)
expect_equal(# same number of best designs
  dim(uniqueBest_gen2)[1],
  length(optSD_gen2$SD)
)
for (i in 1:dim(uniqueBest_gen2)[1]){
  expect_equal(# same best designs
    sort(c(locAll_l[uniqueBest_gen2[i,] == 1], locKeep_l)),
    sort(optSD_gen2$SD[[i]])
  )  
}

#-------- call via optimiseSD -------------------- #
optimiseSD_genetic_m2 = replaceDefault(
    optimiseSD_genetic, newDefaults = list(
      evalFunc = numberPenalty, # should have another function to test if it is forwarded properly
      popSize = 20, 
      iters = 20
    ),
    type = "optimisationFun.optimiseSD")[[1]]

set.seed(07021916)
OptSD_gen1 = optimiseSD(
  simulations = radioactivePlumes_local, 
  costFun = singleDetection,
  optimisationFun = optimiseSD_genetic_m2,
  locationsAll = locAll_l, 
  locationsFix = locKeep_l, 
  locationsInitial = locInit_l[1,],
  aimCost = 0.05,
  nameSave = NA # needed, as else NULL is forwarded and parameter destroyed 
)

optimiseSD_genetic_1 = replaceDefault(
  optimiseSD_genetic, newDefaults = list(
    evalFunc = numberPenalty, 
    popSize = 20, 
    iters = 10
  ),
  type = "optimisationFun.optimiseSD")[[1]]

set.seed(09112006)
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
expect_equivalent(
  OptSD_gen1,
  optSD_gen1
)

# other simulations, costFun...
data(SimulationsSmall)
meanFun = function(x){mean(x, na.rm = TRUE)}
minDist = replaceDefault(
  spatialSpread, newDefaults = list(
    fun = minimalDistance,
    fun_R = meanFun
  ), type = "costFun.optimiseSD"
)[["fun"]] 

optSD_gen_SS = optimiseSD_genetic(
  simulations = SimulationsSmall,
  costFun = minDist,
  aimNumber = 2
  )


#29061929
#18031931
#14012005
#21121993