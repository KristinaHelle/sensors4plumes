################################################################################
#  Spatial Simulated Annealing based on vanGroenigen and intamapInteractive    #
################################################################################
#  - - - - - - - - - - - - help functions - - - - - - - - - - - - - - - - - - - #
# _ - _ - _ - _ - _ - _ - adjustShift - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ #
# adjust shift size parameters according to size of area
adjustShift = function(simulations, locationsAll, minShiftFactor, maxShiftFactor){
  allLocations = subsetSDF(simulations@locations, locations = locationsAll)
#  coordinatesAll = coordinates(allLocations)
  extentAll = apply(X = bbox(allLocations), FUN = diff, MARGIN = 1)
  shift = extentAll %*% t(c(minShiftFactor, maxShiftFactor))
  dimnames(shift)[[1]] = c("x", "y")
  dimnames(shift)[[2]] = c("min", "max")
  # correct if smaller than grid size
  if (is.element(class(simulations@locations),
                 c("SpatialPixelsDataFrame", "SpatialGridDataFrame", "SpatialPolygridDataFrame"))){
    shiftTooSmall = simulations@locations@grid@cellsize > shift[,"min"]
    shift[, "min"][shiftTooSmall] = simulations@locations@grid@cellsize[shiftTooSmall]
    shift[, "max"] = apply(FUN = max, X = shift, MARGIN = 1)
    if (any(shiftTooSmall)){
      warning("Given 'minShiftFactor' is too small, minimal shift is set to cellsize of grid.")
    }
  }
  result = list()
  result[["shift"]] = shift
  result[["allLocations"]] = allLocations
  return(result)
}

# _ - _ - _ - _ - _ - _ - _ - _ - estimateCostDiff - _ - _ - _ - _ - _ - _ - _ - _ - _ #
estimateCostDiff = function(
  simulations,
  locationsInitial,
  locationsFix,
  costFun,
  costInitial,
  n = 5
){
  nAll = nLocations(simulations)
  # sample
  changes = matrix(nrow = 2 * n, ncol = 3)
  changes[,1] = sample.int(length(locationsInitial), 2 * n, replace = TRUE)
  changes[1:n,2] = sample.int(nAll, n, replace = TRUE)
  firstHalf = locationsInitial[changes[n + 1:n,1]] <= nAll/2
  changes[n + 1:n,2][firstHalf] = locationsInitial[changes[n + 1:n,1][firstHalf]] + 1
  changes[n + 1:n,2][!firstHalf] = locationsInitial[changes[n + 1:n,1][!firstHalf]] - 1

  # estimate cost difference
  for (i in 1:(2 * n)){
    changes[i,3] = abs(costInitial - costFun(
      simulations = simulations,
      locations = c(locationsInitial[-changes[i,1]], locationsFix, (1:nAll)[changes[i,2]]))[[1]])
  }
  costDiffMedian = c(median(changes[1:n, 3]), median(changes[n + 1:n, 3]))
  return(costDiffMedian)
}


# _ - _ - _ - _ - _ - _ - moveSensors - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ #
moveSensors = function(
  allLocations,
  locationsPrev,
  maxShiftNumber,
  maxIterations,
  k,
  startMoveProb,
  shift
){
  end = FALSE
  # shift size for this iteration k (decreases linearly)
  shiftSize = shift[,"max"] - apply(FUN = diff, X = shift, MARGIN = 1) * k / maxIterations
  # all coordinates
  coordinatesAll = coordinates(allLocations)

  # which sensors change; which_shifts[1] is always changed;
  #  the others with probability start_q*(1 - k/maxIterations)
  which_shifts = sample.int(length(locationsPrev), size = maxShiftNumber)
  if_shifts = runif(maxShiftNumber) < startMoveProb * (1 - k / maxIterations)
  if_shifts[1] = TRUE

  #locationChange = integer(maxShiftNumber) # may be used for documentation and plotting
  #allNeighbours = list()
  for (i in 1:maxShiftNumber) {
    if (if_shifts[i]){
      coordinatesOld = coordinatesAll[allLocations@data$indexOrig == locationsPrev[which_shifts[i]],]
      #locationsChange = locationsPrev[which_shifts[i]]
      neighbours = abs(coordinatesAll[,1] - coordinatesOld[1]) <= shiftSize["x"] &
        abs(coordinatesAll[,2] - coordinatesOld[2]) <= shiftSize["y"]   # shift within neighbourhood of size x_shift x y_shift; which uniform probability
      neighboursNew = setdiff(allLocations@data$indexOrig[neighbours], locationsPrev)
      if(length(neighboursNew) == 0 | length(coordinatesOld) == 0){# second criterion stops if invalid locationsPrev (should never be needed)
        warning(paste("Terminated after", k, "iterations because allowed shift too small to move sensors."))
        end = TRUE
        break
      }else{
        if (length(neighboursNew) == 1){
          locationsPrev[which_shifts[i]] = neighboursNew # using 'sample' here is 'sample.int' and generates invalid values
        } else {
          locationsPrev[which_shifts[i]] = sample(neighboursNew, size = 1)
        }
      }
#      allNeighbours[[i]] = neighboursNew
    }
    if(end) break
  }
  out = list(
    locationsPrev = locationsPrev,
    end = end
)
  return(out)
}

# -------------------- optimiseSD_ssa ------------------------------------------ #
# cost function is arbitrary but parameters were chosen for plume detection
optimiseSD_ssa = function(
                  simulations,
                  costFun,
                  locationsAll =1:nLocations(simulations),        # all possible sensor locations (as indices of simulations@locations)
                  locationsFix = integer(0),                    # these sensors are always included to determine cost, they cannot be deleted
                  locationsInitial = integer(0),                # current sensors that can be deleted
                  aimCost = NA,                                  # desired cost, stop optimisation if reached
                  aimNumber = NA,                                # not used but needed as parameter if used in optimiseSD
                  nameSave = NA,                                 # without suffix .Rdata
                  plot = FALSE,                                    # if each iteration is plotted
                  verbatim = FALSE,                                 # print intermediate results and keep intermediate steps
                  maxShiftNumber = length(locationsInitial),       # how many of the sensors can be shifted in each iteration
                  startMoveProb = 0.5,                             # initial probability to move (more than one) sensors
                  maxShiftFactor = 0.2,                            # initial max size of shift (as fraction of extent)
                  minShiftFactor = 0.01,                           # final max size of shift
                  maxIterations = 5000,                            # maximal number of iterations (iteration means tried sensor changes, not all are accepted)
                  maxStableIterations = 500,                       # number of sets of sensors to be rejected in a row before it stops
                  maxIterationsJumpBack = maxStableIterations - 1, # if for maxIterationsJumpBack steps it has not become better than ever before, it jumps back to the best result up to now
                  acceptanceMethod = "vanGroenigen",               # = c("intamapInteractive", "vanGroenigen") formula to compute probability to accept worse result; "vanGroenigen" accepts relative to worsening;
                  start_acc_vG = 0.5,                              # acceptance worse locations in first iteration ("vanGroenigen")
                  end_acc_vG = 0.0001,                             # acceptance of worse locations in maxIterations-th iteration ("vanGroenigen")
                  cooling_vG = 0.999,                              # cooling parameter c of vanGroenigen method
                  startAcceptance_vG = 0.5,                        # cooling parameter a1 of vanGroenigen method
                  start_acc_iI = 0.2                               # calibration factor of acceptance for "intamapInteractive" (equals "start_p" in ssaOptim)
  ){
  # parameter check
  if (is.na(aimCost)){
    stop("'optimiseSD_ssa' needs a value for 'aimCost'.")
  }
  if (!is.na(aimNumber)){
    if (aimNumber <= length(locationsFix)){
      stop(paste0("'aimNumber' (", aimNumber, ") must be higher than the number of fix sensors (", length(locationsFix), ")."))
    }
  }
  ###------------- global settings -----------------------------------

  ## acceptance: when to replace old and best design (and reset improvement counter)
  # depends on rank of costCurrent, costOld , costBest;
  # only consider cases with rankGlobal <= rankOld (other ones cannot occur)
  # only consider possible rank patterns
  # rules:
  # Current <= Old:  replace Old design
  # Current < Old:   reset "no improvement"-counter
  # Current <= Best: replace Best design
  # Current < Best:  reset "no improvement of Best"- counter
  ssaTable = data.frame(
    rank          = c(    3,     3,     2,     2,     1,     1,     1,     1),
    rankOld       = c(    1,     2,     2,     3,     1,     3,     2,     3),
    rankGlobal    = c(    1,     1,     1,     1,     1,     1,     2,     2),
    acceptOld     = c(FALSE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE),
    acceptGlobal  = c(FALSE, FALSE, FALSE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE),
    improveOld    = c(FALSE, FALSE, FALSE,  TRUE, FALSE,  TRUE,  TRUE,  TRUE),
    improveGlobal = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,  TRUE,  TRUE)
  )

  ## if not 'locationsInitial' start from random locations
  if (all(is.na(locationsInitial))){
    if (is.na(aimNumber)){
      stop("'locationsInitial' or 'aimNumber' needed.")
    } else {
      locationsInitial = sample(setdiff(locationsAll, locationsFix), aimNumber - length(locationsFix))
      warning(paste0("Starting from ", aimNumber , "(= aimNumber) random locations."))
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

  ## adjust maximal number of shifts per iteration
  maxShiftNumber = min(maxShiftNumber, length(locationsInitial))

  ## initial cost
  costInitial = costFun(simulations = simulations,
                        locations = c(locationsInitial, locationsFix))[[1]]#, plot = FALSE)[[1]]

  ## compute (and adjust) shift size
  shiftSize = adjustShift(simulations, locationsAll, minShiftFactor, maxShiftFactor)
  allLocations = shiftSize[["allLocations"]]
  allLocations@data$indexOrig = locationsAll

  ## determine cooling parameters
  if (acceptanceMethod == "vanGroenigen") {
    costDiffEstimate = estimateCostDiff(simulations = simulations,
                                        costFun = costFun,
                                        locationsInitial = locationsInitial,
                                        locationsFix = locationsFix,
                                        costInitial = costInitial)
    if (all(costDiffEstimate > 0) & !is.na(start_acc_vG) & !is.na(end_acc_vG)){
      cooling = ((log(start_acc_vG) * costDiffEstimate[2]) /
                   (log(end_acc_vG) * costDiffEstimate[1])) ^ {1 / (maxIterations - 1)}
      startAcceptance = - costDiffEstimate[1]/
        (cooling * log(start_acc_vG))
    } else {
      warning("Estimated cost difference when changing one sensor is 0,
              therefore cooling cannot be adapted to fulfil initial and final acceptance rate,
              given values for 'cooling' and 'startAcceptance' are used.")
      cooling = cooling_vG
      startAcceptance = startAcceptance_vG
    }
    cooling = min(cooling, 1)
    startAcceptance = min(startAcceptance, 1)
    # print(paste0("cooling: ", cooling, "; startAcceptance: ", startAcceptance))
  }

  ## monitoring
  report = data.frame(cost = numeric(maxIterations),
                      accepted = logical(maxIterations))
  if (verbatim){
    report = data.frame(costTest = numeric(maxIterations),
                        cost = numeric(maxIterations),
                        costBest = numeric(maxIterations),
                        chi = numeric(maxIterations),
                        accepted = logical(maxIterations))
    SDs = matrix(nrow = maxIterations, ncol = length(locationsInitial) + length(locationsFix))
    SDs_test = matrix(nrow = maxIterations, ncol = length(locationsInitial) + length(locationsFix))
    SDs_best = matrix(nrow = maxIterations, ncol = length(locationsInitial) + length(locationsFix))

    jumpBack = integer(0)
  }
  if (!is.na(nameSave)){
    filename = paste0(nameSave, ".Rdata")
    message(paste0("Sampling designs saved at ", filename, "."))
  }
  ### --------------------------------- initialise ------------------------------
  ## end
  end = FALSE

  ## sensors
  locationsCurrent = c(locationsInitial, locationsFix)#union(locationsInitial, locationsFix)

  locationsOld = locationsCurrent                         # always locations of last accepted iteration
  locationsBest = locationsCurrent                        # globally best locations up to now

  ## cost
  costOld = costInitial                                   # always cost of last accepted iteration
  costBest = costInitial                                  # best cost up to now


  ## iterations
  count = 0                                               # steps since last improvement of costCurrent
  countGlobal = 0                                         # steps since last improvement of global optimum: costBest



  ### ----------------------------- optimisation --------------------------------------
  for (k in 1:maxIterations){
    movedSensors = moveSensors(
      allLocations = allLocations,
      locationsPrev = locationsOld[1:length(locationsInitial)],
      maxShiftNumber = maxShiftNumber,
      maxIterations = maxIterations,
      k = k,
      startMoveProb = startMoveProb,
      shift = shiftSize[["shift"]]
    )
    locationsCurrent = c(movedSensors[["locationsPrev"]], locationsFix)
    end = movedSensors[["end"]]
    if(end) break

    # compute cost of changed sensors
    costCurrent = costFun(simulations = simulations, locations = locationsCurrent)[[1]]#, plot = FALSE)[[1]]
#    report$cost[k] = costOld
#    print(paste0(k, " cost: ", costCurrent[[1]]))
    # decide about acceptance of change
    ranks = rank(c(costCurrent, costOld, costBest), ties.method = "min")
    rankCase = which(apply(FUN = all, X = apply(FUN = "==", X = ssaTable[,1:3], y = ranks, MARGIN = 1), MARGIN = 2))
#    if (verbatim) {
#      report$costTest[k] = costCurrent
#      report$costBest[k] = costBest
#    }

    # if at least as good as last iteration, use new locationsS
    if(ssaTable$acceptOld[rankCase]){ # replace locationsOld
      locationsOld = locationsCurrent
      costOld = costCurrent
      report$accepted[k] = TRUE
    }
    if(!ssaTable$improveOld[rankCase]){ # if no improvement, increase counter
      count = count + 1
    }else{
      count = 0
    }

    # if at least as good as global optimum, replace it
    if(ssaTable$acceptGlobal[rankCase]){ # replace locationsBest
      locationsBest = locationsCurrent
      costBest = costCurrent
    }
    if(!ssaTable$improveGlobal[rankCase]){ # if no improvement, increase counter
      countGlobal = countGlobal + 1
    }else{
      countGlobal = 0
      if (verbatim) print(paste(k, ": global best cost improved to ", signif(costBest, digits = 5), sep = ""))
    }

    # if worse; accept with probability chi
    if(!ssaTable$acceptOld[rankCase]){
      chi = switch(acceptanceMethod,
                   "vanGroenigen" = exp((costOld - costCurrent) / (startAcceptance * cooling ^ k)),
                   "intamapInteractive" = start_acc_iI * exp(-10 * k / maxIterations)
      )
      if (verbatim){
        report$chi[k] = chi
      }
      if(chi > runif(1)){# if accept
        if (verbatim) print(paste0(k, ": worsening of to cost to ", signif(costCurrent, digits = 5), " accepted, chi = ", round(chi, digits = 5), ", costDiff = ", signif(costOld - costCurrent, digits = 5)))
        locationsOld = locationsCurrent
        costOld = costCurrent
        report$accepted[k] = TRUE
      }
    }

    # jump back if no global improvement for maxIterationsJumpBack iterations
    if(countGlobal >= maxIterationsJumpBack){
      locationsOld = locationsBest
      costOld = costBest
      countGlobal = 0

      if (verbatim){
        jumpBack = c(jumpBack, k)
        print(paste(k, ": jump back to optimum of cost ", signif(costOld, digits = 5),  sep = ""))
      }
    }
    # document
    report$cost[k] = costOld
    if (verbatim) {
      SDs_test[k,] = locationsCurrent
      SDs[k,] = locationsOld
      SDs_best[k,] = locationsBest
      report$costTest[k] = costCurrent
      report$costBest[k] = costBest
    }
    # print(paste0(k, " cost: ", report$cost[k]))
    # end if no local improvement for maxStableIterations iterations
    if(count >= maxStableIterations){
      warning(paste("Terminated after", k, "iterations because no improvement for ", maxStableIterations, " iterations.", sep = ""))
      break
    }
  if (!is.na(nameSave) & verbatim){
#    if(k > maxIterations * 0.9){
      save(SDs, report, file = filename)
#    }
  }

    # if desired plot each new, accepted sensor set
#     if (plot & report$accepted[k]) {
#       plotSD(
#         simulations = simulations,
#         costFun = costFun,
#         type = "ssa",
#         SD = locationsOld,
#         param = list(
#           locationsFix = locationsFix,
#           locationsInitial = locationsInitial,
#           locationsChange = locationsChange,
#           locationsPotential = allNeighbours,
#           iterations = k,
#           costBest = costBest)
#         )

#       CostCurrent = costFun(simulations = simulations, locations = locationsOld)
#       if (!is.na(CostCurrent[["image"]])) {
#            require(sp)
#            PlotPtsFix = list("sp.points", coordinates(simulations@locations)[locationsFix, , drop = FALSE],       pch = 2, col = "grey")
#            PlotPtsInit = list("sp.points", coordinates(simulations@locations)[locationsInitial, , drop = FALSE],  pch = 1, col = "grey")
#            PlotPts =  list("sp.points", coordinates(simulations@locations)[locationsOld, , drop = FALSE],     pch = 20, col = "green")
#            PlotPtsPotential = list("sp.points", coordinates(simulations@locations)[allNeighbours, , drop = FALSE],pch = ".", col = "red")
#            PlotPtsChange = list("sp.points", coordinates(simulations@locations)[locationsChange, , drop = FALSE], pch = 4, col = "red")
#            Image = extractSpatialDataFrame(simulations, layer = 1, plume = 1)
#            Image@data$cost = CostCurrent[["image"]]
#            print(spplot(Image, "cost", col.regions = bpy.colors(100), sp.layout = list(PlotPtsInit, PlotPtsFix, PlotPts, PlotPtsPotential,PlotPtsChange),
#                  main = paste("cost = ", signif(costOld, digits = 5),
#                               "\n", "best cost = ", signif(costBest, digits = 5), "\n", "Iterations = ", k)))
#           }
#     }
    if (!is.na(aimCost)) {
      if (costBest <= aimCost) {
        print(paste(k, ": optimisation reached desired cost ", signif(costBest, digits = 5), sep = ""))
        end = TRUE
      }
    }
    if (end) break
  } # end k

  # output
  out = list()
  out[["SD"]] = list()
  out[["SD"]][[as.character(length(locationsBest))]] = matrix(locationsBest, nrow = 1)
  out[["evaluation"]] = data.frame(cost = costBest,
                             number = length(locationsBest))
  out[["report"]] = list()
  out[["report"]][["SD_final"]] = locationsOld
  out[["report"]][["SD_best"]] = locationsBest
  out[["report"]][["cost_final"]] = costOld
  out[["report"]][["cost_best"]] = costBest
  out[["report"]][["iterations"]] = report[1:k,]

  if (verbatim) {
    out[["report"]][["SDs"]] = SDs[1:k,]
    out[["report"]][["SDs_test"]] = SDs_test[1:k,]
    out[["report"]][["SDs_best"]] = SDs_best[1:k,]
    if (acceptanceMethod =="vanGroenigen") {
#      out[["costDiffEstimate"]] = costDiffEstimateAll
      out[["report"]][["coolingPar"]] = c(cooling, startAcceptance)
    }
    out[["report"]][["jumpBack"]] = jumpBack
  }
  return(out)
}
