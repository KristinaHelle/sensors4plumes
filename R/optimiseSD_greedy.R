################################################################################
#                          greedy optimisation                                 #
################################################################################
############ help functions ###################################
# ----------------- find least important sensor to delete ---------------------#
# find the current sensor where sampling design without this sensor is of minimal cost
deleteSensor = function(
  simulations,
  costFun,
  locationsDeletable,      # indices of current sensors that are not fix; one of them shall be deleted
  locationsKeep                 # indices of fix sensors; they are use to compute cost, none of them is deleted
){
  nChangeSensors = length(locationsDeletable)
  costWithoutThisSensor = rep(NA, times = nChangeSensors)
  for(i in 1:nChangeSensors){                                                   # for each current sensor that could be deleted
    locIndicesWithout = c(locationsDeletable[-i], locationsKeep)                     # sampling design without this sensor
    #costFunLoc = replaceDefault(costFun, newDefaults = list(locations = locIndicesWithout))[[1]]
    #costWithoutThisSensor[i] = do.call(what = costFunLoc, args = as.list(formals(costFunLoc)))[["cost"]]
    costWithoutThisSensor[i] = costFun(simulations = simulations, locations = locIndicesWithout)[[1]]#, plot = FALSE)[[1]]      # cost for sampling design without this sensor
  }
  leastImportant = which(costWithoutThisSensor == min(costWithoutThisSensor))   # sensors where cost is minimal without this sensor
  deletableNew = locationsDeletable[-leastImportant[1]]                           # sampling design where this sensor is deleted (delete the first one if there are several)
  return(list(deletableNew, leastImportant, costWithoutThisSensor))
}

# ----------------------- find best sensor to add -----------------------------#
# find the currently empty grid cell where adding a sensor reduces cost most
addSensor = function(
  costFun,
  simulations,
  locationsAddable,        # currently empty potential sensor locations; one of them shall be added
  locationsCurrent         # current sampling design
){
  nAddable = length(locationsAddable)
  costAddedThisSensor = rep(NA, times = nAddable)
  for(i in 1:nAddable){                                                         # for each grid cell where a sensor could be added
    locIndicesAdded = c(locationsAddable[i], locationsCurrent)                          # current sampling design and this sensor added
    #costFunLoc = replaceDefault(costFun, newDefaults = list(locations = locIndicesAdded))[[1]]
    #costAddedThisSensor[i] = do.call(what = costFunLoc, args = as.list(formals(costFunLoc)))[["cost"]]
    costAddedThisSensor[i] = costFun(simulations = simulations, locations = locIndicesAdded)[[1]]#, plot = FALSE)[[1]]          # cost for this sampling design
  }
  best = which(costAddedThisSensor == min(costAddedThisSensor))                 # sensors where cost is minimal when they (i.e. one of them) are added
  Add = locationsAddable[best[1]]                                                 # this best sensor (the first one if there are several)
  return(list(Add, best, costAddedThisSensor))
}

#---------------- determine if next step is add or delete ---------------------#
# determine if sensors should be added or deleted based on comparing current and desired cost / number of sensors
determineNextStep = function(
  locationsCurrent,                      # indices of current sensors, including fix ones
  locationsFix,                      # indices of fix sensors
  locationsAll,            # indices of all possible sensor locations, including fix ones
  kindAim,                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
  aim,                                        # desired number of sensors or desired cost
  maxIterations,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria
  l,                                            # iteration step
  costFun,                                # to be forwarded to computeCost
  simulations
){
  #costFunLoc = replaceDefault(costFun, newDefaults = list(locations = locationsCurrent))[[1]]
  #costHere = do.call(what = costFunLoc, args = as.list(formals(costFunLoc)))[["cost"]]
  costHere = costFun(simulations = simulations, locations = locationsCurrent)[[1]]#, plot = FALSE)[[1]]
  numberHere = length(locationsCurrent)
  switch(kindAim,
         "number" = {
           aimDifference = (numberHere - aim)   # Are there less sensors than desired? (if aim is certain number of sensors)
         },
         "cost" = {
           aimDifference = (aim - costHere[[1]])     # Is cost too high? (if aim is cost below a certain threashold)
         }
  )
  if(aimDifference < 0){                             # => then add sensor; else delete
    NextStep = "add"
  }
  if (aimDifference > 0){
    NextStep = "del"
  }
  if (aimDifference == 0){
    NextStep = "finish"
  } else {
    # exceptions: return message and go to "stop"
    if(l > maxIterations){                                                                    # if maximal number of iterations is exceeded
      NextStep = "stop"
      warning(paste("Optimisation stopped after ", maxIterations, " iterations without reaching a stable result!", sep = ""))
    }
    if (length(setdiff(locationsCurrent, locationsFix)) == 0 & aimDifference > 0){# this case should never occur as locationsFix should be part of locationsCurrent
      NextStep = "stop"
      switch(kindAim,
             "cost" = { warning("All but the fix sensors are deleted and cost is still below 'aimCost'.")  },
             "number" = {warning("All but the fix sensors are deleted but number of sensors is still above 'aimNumber'.")}
      )
    }
    if (length(setdiff(locationsAll, locationsCurrent)) == 0 & aimDifference < 0){
      NextStep = "stop"
      switch(kindAim,
             "cost" = { warning("All possible sensors are included but cost is still above 'aimCost'.")  },
             "number" = {warning("All possible sensors are included and number of sensors is still below 'aimNumber'.")}
      )
    }
  }
#   if(costHere[[1]] == 0 & aimDifference < 0 & kindAim == "number"){                       # if cost is 0 with less than the desired number of sensors
#     NextStep = "stop"
#     warning(paste("Already ", numberHere ," sensors are sufficient to get cost = 0. The Optimisation was stopped at the first sampling design with cost = 0, there may be such sampling designs with less sensors.", sep = ""))
#   }
#
#   if(length(setdiff(locationsCurrent, locationsFix)) == 0  & aimDifference >= 0 & kindAim == "cost"){   # if no more sensors can be deleted and the desired cost is not exceeded
#     NextStep = "stop"
#     warning("The desired cost can be reached without any sensors / with only the fixed sensors.")
#   }
#
#   if(length(setdiff(locationsAll, locationsCurrent)) == 0 & aimDifference < 0 & kindAim == "cost"){# if no more sensors can be added but desired cost is exceeded (this should only happen if aim is negative what can never be achieved)
#     NextStep = "stop"
#     warning("The desired cost cannot be achieved, even with sensors at all 'locationsAll'.")
#   }
#
#   if(length(setdiff(locationsCurrent, locationsFix)) == 0  & aimDifference > 0 & kindAim == "number"){   # if no more sensors can be deleted and the desired cost is not exceeded
#     NextStep = "stop"
#     warning("The desired number is below the number of fix sensors.")
#   }

  # message
  switch(kindAim,
         "number" = {
           message(paste("iteration: ", l, "; number of sensors: ", numberHere, " (aim: ", aim, "); cost: ", costHere[[1]], "; next step: ", NextStep))
         },
         "cost" = {
           message(paste("iteration: ", l, "; number of sensors: ", numberHere, "; cost: ", costHere[[1]],  " (aim: ", aim, "); next step: ", NextStep))
         }
  )
  nextStep = list(
    nextStep = NextStep,
    cost = costHere[[1]])
  return(nextStep)
}



optimiseSD_greedy = function(
  simulations,                     # simulations object to optimise sensors for; forwarded to 'costFun', 'locations...' must be subset of locations else they are ignored
  costFun,                         # function to compute cost
  locationsAll =1:nLocations(simulations),        # all possible sensor locations (as indices of simulations@locations)
  locationsFix = integer(0),                    # these sensors are always included to determine cost, they cannot be deleted
  locationsInitial = integer(0),                # current sensors that can be deleted
  aimCost = NA,
  aimNumber = NA,
  nameSave = NA,
  plot = FALSE,
  verbatim = FALSE,
  #  kindAim,                         # c("cost", "number"): if aim is "cost", sensors are added until cost is below 'aim'; if aim is "number", sensors are reduced to not extend this number
  #  aim,                             # list: number, cost; cost limit or maximal number of sensors;
  #  keepInitial = FALSE,                     # if initial sensors may be replaced
  maxIterations = 100,             # after this number of iterations search is stopped
  swap = FALSE                      # if optimisation is to be continued by adding and deleting in turns after greedy optimum is reached
){
  costFun = replaceDefault(costFun, newDefaults = list(simulations = simulations))[[1]]
  # ------------------ locations -------------------------
  locations = cleanIndices(
    locationsTotal = 1:nLocations(simulations),
    locationsAll = locationsAll,
    locationsFix = locationsFix,
    locationsInitial = locationsInitial
  )
  locationsAll = locations[["locationsAll"]]
  locationsFix = locations[["locationsFix"]]
  locationsInitial = locations[["locationsInitial"]]
  locationsCurrent = union(locationsInitial, locationsFix)

  # ------------------- aim --------------------------------
  if (!is.na(aimNumber)){
    if (aimNumber <= length(locationsFix)){
      stop(paste0("'aimNumber' (", aimNumber, ") must be higher than the number of fix sensors (", length(locationsFix), ")."))
    }
  }

  whichAim = c(!is.na(aimNumber), !is.na(aimCost))
  if (identical(whichAim, c(TRUE, TRUE))){
    warning("The optimisation will aim to achieve the desired cost 'aimCost';
            the desired number of sensors 'aimNumber' is ignored.")
    kindAim = "cost"
    aim = aimCost
  }
  if (identical(whichAim, c(TRUE, FALSE))){
    kindAim = "number"
    aim = aimNumber
  }
  if (identical(whichAim, c(FALSE, TRUE))){
    kindAim = "cost"
    aim = aimCost
  }
  if (identical(whichAim, c(FALSE, FALSE))){
    warning("No aim for the optimisation indicated; for the given number of 'locationsInitial' cost is minimised.")
    kindAim = "number"
    aim = length(locationsInitial) + length(locationsFix)
  }

  if (!is.na(nameSave)){
    filename = paste0(nameSave, ".Rdata")
    message(paste0("Sampling designs saved at ", filename, "."))
  }

  ############ main algorithm ###################################
  # initialise

  SDs = list()
  l = 1
  SDsDiffer = TRUE
  SDs[[1]] = locationsCurrent

  # compare to aim, decide if sensors must be added or deleted
  NEXT = determineNextStep(
    locationsCurrent = locationsCurrent,
    locationsFix = locationsFix,
    locationsAll = locationsAll,
    kindAim = kindAim,
    aim = aim,
    maxIterations = maxIterations,
    l = l,
    costFun = costFun,
    simulations = simulations
  )

  if(NEXT[[1]] == "del"){ # delete sensors because there are too many / cost is below the aim
    #print("loop primary del")
    # initialise from which sets sensors shall be added or deleted
    locationsKeep = locationsFix
    AddFrom = setdiff(locationsAll, locationsKeep)

    message("As there are more sensors than necessary, sensors are deleted one by one.
            Always deleting the least important of the (initally) ", length(setdiff(locationsCurrent, locationsKeep)), " non-fix sensors.")

    while(SDsDiffer & NEXT[[1]] != "stop"){ # stop, if same SD it repeated
      #print(paste("loop del SDsDiffier AND no stop",SDsDiffer, NEXT[[1]]))
      # delete sensors
      locationsDeletable = setdiff(locationsCurrent, locationsKeep)

      while(NEXT[[1]] == "del"){
        #print("loop inner del")
        # determine which sensor to delete
        DeletableNew = deleteSensor(
          simulations = simulations,
          costFun = costFun,
          locationsDeletable = locationsDeletable,
          locationsKeep = locationsKeep
        )
        # delete sensor: update sensors that can be deleted and current ones
        locationsDeletable = DeletableNew[[1]]
        locationsCurrent = c(locationsDeletable, locationsKeep)
        l = l + 1
        # save current sampling design
        SDs[[l]] = locationsCurrent
        if (!is.na(nameSave)){
          save(SDs, file = filename)
        }


        # is aim reached?
        NEXT = determineNextStep(
          locationsCurrent = locationsCurrent,
          locationsFix = locationsFix,
          locationsAll = locationsAll,
          kindAim = kindAim,
          aim = aim,
          maxIterations = maxIterations,
          l = l,
          costFun = costFun,
          simulations = simulations
        )
        if (NEXT[[1]] == "add" | NEXT[[1]] == "finish"){
          message("Greedy optimum found.")

          if (swap){
            message("As 'swap = TRUE' optimisation is continued by adding and deleting sensors in turns.")
          }
      }
    }
    if (swap & NEXT[[1]] != "stop"){# continue by adding sensors (and then deleting them in turns)
      #print("if swap")
      # add sensors
      locationsAddable = setdiff(AddFrom, locationsCurrent)
      NEXT[[1]] = "add"
      while(NEXT[[1]] == "add"){
        #print("loop inner add")
        # decide which sensor to add
        Add = addSensor(
          costFun = costFun,
          simulations = simulations,
          locationsAddable = locationsAddable,
          locationsCurrent = locationsCurrent
        )
        # current sampling design
        locationsCurrent = c(Add[[1]], locationsCurrent)
        locationsAddable = setdiff(AddFrom, locationsCurrent)
        l = l + 1
        SDs[[l]] = locationsCurrent
        if (!is.na(nameSave)){
          save(SDs, file = filename)
        }


        # determine how to continue: add, delete, or stop
        NEXT = determineNextStep(
          locationsCurrent = locationsCurrent,
          locationsFix = locationsFix,
          locationsAll = locationsAll,
          kindAim = kindAim,
          aim = aim,
          maxIterations = maxIterations,
          l = l,
          costFun = costFun,
          simulations = simulations
        )
      }
    }
    # is sampling design repeated?
    SDsDifferVector = sapply(FUN = setequal, X = SDs, SDs[[l]])
    SDsDiffer = !any(SDsDifferVector[-l])
    if (!swap){
      SDsDiffer = FALSE
    }
    }
    #print(paste("loop end, SDsDiffer: ", SDsDiffer))
  }
  else{ # add sensors because there are too few / cost is above threshold
    #print("loop primary add")
    # initialise from which sets sensors shall be added or deleted
    locationsKeep = locationsFix
    AddFrom = setdiff(locationsAll, locationsKeep)

    message(paste("As there are less sensors than necessary, sensors are added one by one.
                  Always adding one at the most important of the (initally) ", length(setdiff(locationsAll, locationsCurrent)), " empty locations."))

    while(SDsDiffer & NEXT[[1]] != "stop"){ # stop, if same SD it repeated
      #print(paste("loop add SDsDiffier AND no stop",SDsDiffer, NEXT[[1]]))
      locationsAddable = setdiff(AddFrom, locationsCurrent)

      while(NEXT[[1]] == "add"){         # add sensors
        #print("loop inner add")
        # determine which sensor to add
        Add = addSensor(
          costFun = costFun,
          simulations = simulations,
          locationsAddable = locationsAddable,
          locationsCurrent = locationsCurrent
        )
        # update current sampling design and save it
        locationsCurrent = c(Add[[1]], locationsCurrent)
        locationsAddable = setdiff(AddFrom, locationsCurrent)
        l = l + 1
        SDs[[l]] = locationsCurrent



        # is aim reached?
        NEXT =  determineNextStep(
          locationsCurrent = locationsCurrent,
          locationsFix = locationsFix,
          locationsAll = locationsAll,
          kindAim = kindAim,
          aim = aim,
          maxIterations = maxIterations,
          l = l,
          costFun = costFun,
          simulations = simulations
        )
        if (NEXT[[1]] == "del" | NEXT[[1]] == "finish"){
          message("Greedy optimum found.")

          if (swap){
            message("As 'swap = TRUE' optimisation is continued by deleting and adding sensors in turns.")
          }
      }
    }

    if (swap & NEXT[[1]] != "stop"){
      #print("if swap")
      # delete sensors
      locationsDeletable = setdiff(locationsCurrent, locationsKeep)
      NEXT[[1]] = "del"
      while(NEXT[[1]] == "del"){
        #print("loop inner del")
        # decide which sensor to delete
        DeletableNew = deleteSensor(
          simulations = simulations,
          costFun = costFun,
          locationsDeletable = locationsDeletable,
          locationsKeep = locationsKeep
        )
        # update and save current sampling design
        locationsDeletable = DeletableNew[[1]]
        locationsCurrent = c(locationsDeletable, locationsKeep)
        l = l + 1
        SDs[[l]] = locationsCurrent
        if (!is.na(nameSave)){
          save(SDs, file = filename)
        }
        # determine how to continue: add, delete, or stop
        NEXT =  determineNextStep(
          locationsCurrent = locationsCurrent,
          locationsFix = locationsFix,
          locationsAll = locationsAll,
          kindAim = kindAim,
          aim = aim,
          maxIterations = maxIterations,
          l = l,
          costFun = costFun,
          simulations = simulations
        )
      }
    }
    # is sampling design repeated?
    SDsDifferVector = sapply(FUN = setequal, X = SDs, SDs[[l]])
    SDsDiffer = !any(SDsDifferVector[-l])
    if (!swap){
      SDsDiffer = FALSE
    }
    }
    #print(paste("loop end, SDsDiffer: ", SDsDiffer))
  }

  # summarise result
  costSDs = numeric(length(SDs))
  for (i in seq(along = costSDs)){
    costFunLoc = replaceDefault(costFun, newDefaults = list(locations = SDs[[i]]))[[1]]
    costSDs[i] = do.call(what = costFunLoc, args = as.list(formals(costFunLoc)))[["cost"]]
  }
  #costSDs =
  #  sapply(lapply(FUN = costFun, X = SDs, simulations = simulations), "[[", "cost")
  lengthSDs = sapply(FUN = length, X = SDs)

  nSDs = unique(sort(lengthSDs))
  SD = list()
  evaluation = data.frame(cost = numeric(length(nSDs)), number = nSDs)
  for (i in seq(along = nSDs)){
    whichN_i = which(lengthSDs == nSDs[i])
    best_i = which(costSDs[whichN_i] == min(costSDs[whichN_i]))
    SD_i = matrix(nrow = length(best_i), ncol = nSDs[i])
    for (j in seq(along = best_i)){
      SD_i[j,] = sort(SDs[[whichN_i[best_i][j]]])
    }
    SD[[i]] = unique(SD_i)
    evaluation$cost[i] = min(costSDs[whichN_i])
  }

  # switch(kindAim,
  #        "number" = {
  #          lengthAim = which(lengthSDs <= aim)
  #          minCost_lengthAim = which(costSDs == min(costSDs[lengthAim]))
  #          shortest_minCost_lengthAim = which(lengthSDs  == min(lengthSDs[intersect(lengthAim, minCost_lengthAim)]))
  #          finalSDwhich = intersect(shortest_minCost_lengthAim, intersect(minCost_lengthAim , lengthAim))
  #        },
  #        "cost" = {
  #          costAim = which(costSDs <= aim)
  #          minLength_costAim = which(lengthSDs == min(lengthSDs[costAim]))
  #          cheapest_minLength_costAim = which(costSDs == min(costSDs[minLength_costAim]))
  #          finalSDwhich = intersect(costAim, intersect(minLength_costAim, cheapest_minLength_costAim))
  #        }
  # )
  # finalSDs = matrix(nrow = length(finalSDwhich), ncol = length(SDs[[finalSDwhich[1]]]))
  #
  if (!is.na(nameSave)){
    save(SD, file = filename) # changed, was SDs
  }
  out = list()
  out[["SD"]] = SD
  out[["evaluation"]] = evaluation
  out[["report"]] = list()
  out[["report"]][["SDs"]] = SDs
  out[["report"]][["evalSDs"]] = data.frame(cost = costSDs, number = lengthSDs)
  # out = list()
  # out[["SD"]] = SDs[[finalSDwhich[1]]]
  # out[["SDs"]] =  SDs
  # out[["finalSDwhich"]] = finalSDwhich
  # out[["evalSDs"]] = data.frame(cost = costSDs, number = lengthSDs)
  return(out)
  }# end optimiseSD_greedy
