################################################################################
#                         genetic optimisation                                 #
################################################################################
############ help functions ###################################
connection_logCost = function(cost, number, weight = 1){
  out = log(cost) + exp(weight * number)
}
connection_penalty = function(cost, number, penalty = 2, aimNumberNF = 1){
  if (number > aimNumberNF){
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
  aimNumberNF,
  penalty = 2
  #plot = FALSE
  ){

  locations = c(locationsAll[chromosome == 1], locationsFix)

  # compute cost for this chromosome
  number = sum(chromosome) # fix locations do not count
#  out = list()
  if (number > aimNumberNF){
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

# old
# rbga2SD = function(
#   simulations,
#   costFun,
#   locationsFix,
#   locationsAll,
#   rbgaResults
#   ){
#   # extract all best SDs
# #  bestChromosomesIndex = which(rbgaResults$evaluations == min(rbgaResults$evaluations))
# #  bestChromosomes = rbgaResults$population[bestChromosomesIndex, , drop = FALSE]
#
#   # delete multiples
#   uniqueChromosomes = unique(rbgaResults$population)
#   nSD = nrow(uniqueChromosomes)
#   cost = numeric(nSD)
#   SD = list()
#   for (i in 1:nSD){
#     chromosomesTotal = rep(0, nLocations(simulations))
#     chromosomesTotal[locationsFix] = 1
#     chromosomesTotal[locationsAll] = uniqueChromosomes[i,]
#     SD[[i]] = which(chromosomesTotal == 1)
#     cost[i] = costFun(simulations, locations = SD[[i]])[[1]]
#   }
#   out = list(SD = SD, cost = cost)
#   return(out)
# }
rbga2SD = function(
  simulations,
  costFun,
  locationsFix,
  locationsAll,
  rbgaResults
){

  # delete multiples
  uniqueChromosomes = unique(rbgaResults$population)
  nSD = nrow(uniqueChromosomes)
  # sort size
  numberSD = apply(X = uniqueChromosomes, FUN  = sum, MARGIN = 1)
  sizes = data.frame(table(numberSD))

  eval = data.frame(cost = numeric(nrow(sizes)), number = integer(nrow(sizes)))
  SD = list()
  for (i in seq(along = sizes$Freq)){
    whichChr_i = which(numberSD == sizes$numberSD[i])
    chromosomesTotal_i = matrix(0, nrow = sizes$Freq[i], ncol = nLocations(simulations))
    chromosomesTotal_i[,locationsFix]  = 1
    SDs_i = list()
    cost_i = numeric(sizes$Freq[i])
    for (j in 1:sizes$Freq[i]){
      chromosomesTotal_i[j,locationsAll] = uniqueChromosomes[whichChr_i[j],]
      SDs_i[[j]] = which(chromosomesTotal_i[j,] == 1)
      cost_i[j] = costFun(simulations, locations = SDs_i[[j]])[[1]]
    }
    # extract the best
    best_i = which(cost_i == min(cost_i))

    SD[[i]] = matrix(nrow = length(best_i), ncol = length(SDs_i[[1]]))
    for (j in seq(along = best_i)){
      SD[[i]][j,] = SDs_i[[best_i[j]]]
    }
    eval$cost[i] = min(cost_i)
    eval$number[i] = dim(SD[[i]])[2]
  }
  out = list(SD = SD,
             evaluation = eval,
             report = rbgaResults)
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
  aimNumber = NA,
  nameSave = NA,
  plot = FALSE,
  verbatim = FALSE,
  evalFunc = numberPenalty,
  #monitorFunc = NA,
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
  ### aimNumberFN - without fix sensors
  if (!is.na(aimNumber)){
    aimNumberNF =  aimNumber - length(locationsFix)
    if (aimNumberNF <= 0){
      stop(paste0("'aimNumber' (", aimNumber, ") must be higher than the number of fix sensors (", length(locationsFix), ")."))
    }
  }

  ## generate locationsInitial - if missing
  if (all(is.na(locationsInitial))){
    if (! is.na(aimNumber)){
      locationsInitial = t(replicate(popSize - 1,
                               sample(locationsAll, aimNumberNF)))
      message("All initial locations are chosen randomly.")
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

  if (is.na(aimNumber)){
    if (is.null(dim(locationsInitial))){
      aimNumberNF = length(locationsInitial)
    } else {
      aimNumberNF = length(locationsInitial[1,])
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
      aimNumberNF = aimNumberNF
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
                         zeroToOneRatio = round(nAll/aimNumberNF),
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

  # out = list()
  # out[["SD"]] = list()
  # out[["SD"]][[as.character(length(locationsBest))]] = matrix(locationsBest, nrow = 1)
  # out[["evaluation"]] = data.frame(cost = costBest,
  #                                  number = length(locationsBest))
  # out[["SDfix"]] = locationsFix
  # out[["SDall"]] = locationsAll
  # out[["report"]] = list()
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
