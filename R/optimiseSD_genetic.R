################################################################################
#                         genetic optimisation                                 #
################################################################################
############ help functions ###################################
connection_logCost = function(cost, number, weight = 1){
  out = log(cost) + exp(weight * number)
}
connection_penalty = function(cost, number, penalty = 2, aimNumber = 1){
  if (number > aimNumber){
    out = penalty
  } else {
    out = cost
  }
  return(out)
}

evaluation = function(
  chromosome,
  simulations,
  costFun,
  locationsAll,
  locationsFix,
  connectionFun
  ){
  locations = c(locationsAll[chromosome == 1], locationsFix)
  cost = costFun(simulations = simulations, locations = locations)[["cost"]]
  number = sum(chromosome)  # fix locations do not count
  out = connectionFun(cost = cost,
                      number = number)
  return(out)
}  



numberPenalty = function(
  chromosome,
  simulations,
  costFun,
  locationsAll,
  locationsFix,
  aimNumber,
  penalty = 2,
  plot = FALSE){
  
  locations = c(locationsAll[chromosome == 1], locationsFix)
  
  # compute cost for this chromosome
  number = sum(chromosome) # fix locations do not count
#  out = list()
  if (number > aimNumber){
    #out[["cost"]] = penalty
    out = penalty
  } else {
    cost = costFun(simulations = simulations, locations = locations)
    out = cost[[1]]
    #out[["cost"]] = cost[[1]]
#    if (plot & !is.null(cost[["map"]])){
#       out[["map"]] = cost[["map"]]
#    }
  }
  return(out)
}

rbga2SD = function(
  simulations, 
  costFun,
  locationsFix,
  locationsAll,
  rbgaResults
  ){
  # extract all best SDs 
#  bestChromosomesIndex = which(rbgaResults$evaluations == min(rbgaResults$evaluations)) 
#  bestChromosomes = rbgaResults$population[bestChromosomesIndex, , drop = FALSE]
  
  # delete multiples
  uniqueChromosomes = unique(rbgaResults$population)
  nSD = nrow(uniqueChromosomes)
  cost = numeric(nSD)
  SD = list()
  for (i in 1:nSD){
    chromosomesTotal = rep(0, nLocations(simulations))
    chromosomesTotal[locationsFix] = 1
    chromosomesTotal[locationsAll] = uniqueChromosomes[i,]
    SD[[i]] = which(chromosomesTotal == 1)
    cost[i] = costFun(simulations, locations = SD[[i]])[[1]]
  }
  out = list(SD = SD, cost = cost)
  return(out)
}

# monitor func: user-given?
# could do: saving; stopping when aimCost reached?; plotting?
# these functions could be combined with user-given function body


optimiseSD_genetic = function(
  simulations, 
  costFun, 
  locationsAll = 1:nLocations(simulations), 
  locationsFix = integer(0), 
  locationsInitial = integer(0), # vector or matrix; if matrix: row ~ design
  aimCost = NA, 
  aimNumber = NULL,
  nameSave = NA,
  verbatim = FALSE,
  evalFunc = numberPenalty,
  #monitorFunc = NULL,
  popSize = 200, 
  iters = 100, 
  mutationChance = NA,
  elitism = NA
){
  #require(genalg)
  if (!is.na(nameSave)){
    verbatim = TRUE
  }
  if (!is.na(aimCost)){
    verbatim = TRUE
    if (is.na(nameSave)){
      nameSave = "opt_genetic"
      warning("As 'aimCost' is given, the algorithm will break when this aim is reached without returning a result.
              Therefore all intermediate results are saved at 'opt_genetic.Rdata'")
    }
  }
  ## generate locationsInitial - if missing
  if (missing(locationsInitial)){
    if (! missing(aimNumber)){
      locationsInitial = t(replicate(popSize - 1, 
                               sample(locationsAll, aimNumber))) 
      warning("All initial locations are chosen randomly.")
    } else {
      stop("At least one of 'locationsInitial' or 'aimNumber' must be given.")
    }
  }
 
  ## clean locations
  locations = cleanIndices(
    locationsTotal = 1:nLocations(simulations),
    locationsAll = locationsAll,
    locationsFix = locationsFix,
    locationsInitial = locationsInitial
  ) 
  locationsAll = locations[["locationsAll"]]
  locationsFix = locations[["locationsFix"]]
  locationsInitial = locations[["locationsInitial"]]
  nAll = length(locationsAll)
  nFix = length(locationsFix)

  if (is.null(aimNumber)){
    if (is.null(dim(locationsInitial))){
      aimNumber = length(locationsInitial)  
    } else {
      aimNumber = length(locationsInitial[1,]) 
    }  
  }
  
  cost0sensors = costFun(simulations = simulations, locations = integer(0))

  # prepare evalFunc
  evalFunction = replaceDefault(
    fun = evalFunc,
    newDefaults = list(
      simulations = simulations,
      costFun = costFun,
      locationsAll = locationsAll,
      locationsFix = locationsFix,
      aimNumber = aimNumber
      #plot = plot
    ), type = "evalFunc.rbga.bin"
  )[[1]]
  
  # prepare monitoring
  if (!is.na(nameSave)){
    filename = paste0(nameSave, ".Rdata")
    message(paste0("Populations are saved at ", filename, "."))
    populations = array(dim = c(popSize, nAll, iters))
    save(populations, file = filename)  
  }
  monitor = replaceDefault(monitor,
                              newDefaults = list(
                                simulations = simulations,
                                costFun = costFun,
                                locationsAll = locationsAll,
                                locationsFix = locationsFix,
                                aimCost = aimCost,
                                nameSave = nameSave  
                                ))[[1]]
  
  # chromosome must contain only the locations that can be changed, not fix, not forbidden
  # suggestions, from locationsInitial
  if (is.null(dim(locationsInitial))){
    suggestionsTotal = matrix(0, ncol = length(locationsAll), nrow = 1)
    suggestionsTotal[1, is.element(locationsAll, locationsInitial)] = 1
  } else {
    suggestionsTotal = matrix(0, ncol = length(locationsAll), nrow = dim(locationsInitial)[1])
    for (i in 1:dim(locationsInitial)[1]){
      suggestionsTotal[i, is.element(locationsAll, locationsInitial[i,])] = 1
    }
  }  
  # nrow(suggestions) must not exceed popSize - 1, delete suggestions if necessary
  if (nrow(suggestionsTotal) >= popSize){
    warning(paste0("The number (row) of 'locationsInitial' must be below 'popSize' (=", popSize, "); extra rows are deleted."))
    suggestionsTotal = suggestionsTotal[1:(popSize - 1), , drop = FALSE]
  }

  rbgaResults = rbga.bin(size = nAll, 
                         suggestions = suggestionsTotal,
                         evalFunc = evalFunction,
                         zeroToOneRatio = round(nAll/aimNumber),
                         monitorFunc = monitor,
                         verbose = verbatim,
                         popSize = popSize, 
                         iters = iters, 
                         mutationChance = mutationChance,
                         elitism = elitism)
  
  optimalSD = rbga2SD(simulations = simulations, 
                     locationsFix = locationsFix,
                     locationsAll = locationsAll,
                     rbgaResults = rbgaResults, 
                     costFun = costFun)
  optimalSD[["report"]] = rbgaResults
  # aimNumber: penalty, if too much
  # aimCost: stop if reached
  return(optimalSD)
}
    




monitor = function(results,
                   simulations,
                   costFun,
                   locationsAll,
                   locationsFix,
                   aimCost = NA,
                   nameSave = NA){
  #print(results$evaluations) 
  if (!is.na(nameSave)){
    filename = paste0(nameSave, ".Rdata")
    load(filename)
    populations[,,results$iter] = results$population
    save(populations, file = filename)    
  }
  if (!is.na(aimCost)){
    if(min(results$evaluations) <= aimCost){
      stop(paste0("'aimCost' (", aimCost, ") reached; best cost is ", min(results$evaluations)))      
    }
  }
}

#   monitorFunc = function(rbgaResults){
#     out = plotSD(
#       simulations = simulations,
#       costFun = costFun,
#       type = "genetic", 
#       SD = locationsOld,
#       plot = plot,
#       param = list(
#         locationsFix = locationsFix,
#         locationsInitial = locationsInitial,
#         locationsAll = locationsAll,
#         rbgaResults = rbgaResults
#       )
#     )  
#   }

#
# if (!is.null(nameSave)){
#   monitorFunction = function(rbgaResults, saveName = nameSave){
#     # turn result into SDs and save
#     # if there is a monitor function, also do what it requires
#   }
# } else {
#   if (!is.null(monitorFunc)){
#     monitorFunction = monitorFunc
#   }
# }

## how monitorFunc is called by rgba.bin
## note: if verbose == FALSE, it is not called at all
# if (!is.null(monitorFunc)) {
#   if (verbose) 
#     cat("Sending current state to rgba.monitor()...\n")
#   result = list(type = "binary chromosome", size = size, 
#                 popSize = popSize, iter = iter, iters = iters, 
#                 population = population, elitism = elitism, 
#                 mutationChance = mutationChance, evaluations = evalVals, 
#                 best = bestEvals, mean = meanEvals)
#   class(result) = "rbga"
#   monitorFunc(result)
# }
