# if findSensorNumber = FALSE and returned SD detects all plumes: it is unknown if less sensors can do that as well
# if findSensorNumber = TRUE (and findAllOptimal): the returned number of sensors is minimal to detect all plumes
# if findAllOptima = TRUE: all SDs of given length with maximal detection are returned OR all SDs that can detect all plumes (if they have less sensors)

optimiseSD_global = function(
  simulations,
  costFun = NA,
  locationsAll = 1:nLocations(simulations),
  locationsFix = integer(0),
  locationsInitial = integer(0),
  aimCost = NA,
  aimNumber = NA,
  nameSave = NA,
  plot = FALSE,
  verbatim = FALSE,
  detectable = 1,          # layer of simulations@values with detectability
  maxIterations = NA,         # maximal number of iterations, if no value given: number of all possibilities
  findAllOptima = FALSE,     # if TRUE a second search is run to find all optimal designs
  findSensorNumber = FALSE  # only relevant if findAllOptima = FALSE (else automatically TRUE): if at least 'aimNumber' sensors can detect all plumes, find out how many sensors are needet to detect all plumes?
  ){
  # --------------- prepare other parameters ----------------- #
  locations = cleanIndices(
    locationsTotal = 1:nLocations(simulations),#length(simulations@locations),
    locationsAll = locationsAll,
    locationsFix = locationsFix,
    locationsInitial = locationsInitial
  )
  locationsAll = locations[["locationsAll"]]
  locationsFix = locations[["locationsFix"]]
  locationsInitial = locations[["locationsInitial"]]

  if (!is.na(aimNumber)){
    if (aimNumber <= length(locationsFix)){
      stop(paste0("'aimNumber' (", aimNumber, ") must be higher than the number of fix sensors (", length(locationsFix), ")."))
    }
  }

  if (is.na(aimNumber)){
    aimNumber = length(locationsInitial) + length(locationsFix)
  }
  if (aimNumber <= length(locationsFix)) {
    stop("'aimNumber' (> 0) needed or 'locationsInitial' to use their length as 'aimNumber'.")
  }
  aimNumberNF = aimNumber - length(locationsFix)

  # -------- prepare 'detectable' values ---------- #
  # extract layer of interest
  if (!is.na(nameSave)){
    detectableSim = subset(x = simulations,
                           kinds = detectable,
                           nameSave = paste0(nameSave, "_detectableSim"))
  } else {
    detectableSim = subset(x = simulations,
                           kinds = detectable)
  }

  # check, if layer is logical or consists of 0 and 1 only
  is_01 = FALSE
  if (dataType(detectableSim@values) == "LOG1S"){# data is logical -> TRUE
    is_01 = TRUE
  } else {# check value-wise
    is01 = function(x, nout = 1){
      if (!is.finite(x[1])){# NA, NaN, Inf -> FALSE
        is01 = FALSE
      } else {
        if (x[1] == 1 | x[1] == 0){
          is01 = TRUE
        } else {# numbers except 0 1 -> FALSE
          is01 = FALSE
        }
      }
      return(is01)
    }
    is_detectable = simulationsApply(
      simulations = detectableSim,
      values = 1,
      fun_pl = is01,
      fun_Rpl_cellStats = "sum")
    if (is_detectable$cellStats_global_locationsplumes == ncell(detectableSim@values)){
      is_01 = TRUE
    }
  }
  if (!is_01){
    stop("Global search cannot be applied as the given layer of simulations values is not logical or consisting of 0 and 1.")
  }

  # if locationsFix given, delete all these locations and all plumes detected there
  detectedByFix = rep(0, nPlumes(simulations))
  if (length(locationsFix) > 0){
    # which plumes are detected?
    #if (dataType(detectableSim@values) == "LOG1S"){
    #  detectedByFix = summaryPlumes(simulations = detectableSim,
    #                                locations = locationsFix,
    #                                fun = any)
    #  detectedByFix01 = rep(0, nPlumes(detectableSim))
    #  detectedByFix01[detectedByFix[[2]]] = 1
    #} else {
      detectedByFix = summaryPlumes(simulations = detectableSim,
                                    locations = locationsFix,
                                    values = detectable,
                                    fun = max)[[2]]
    #}
    nDetectedFix = sum(detectedByFix)
    if (nDetectedFix > 0){
      paste0(nDetectedFix, "plume(s) are detected at 'locationsFix' and is/are ignored.")
    }
  }
  detectedByAll = rep(1, nPlumes(simulations))
  if (length(locationsAll) < nLocations(simulations)){
    # which plumes cannot be detected at locationsAll?
    #if (dataType(detectableSim@values) == "LOG1S"){
    #  detectedByAll = summaryPlumes(simulations = detectableSim,
    #                                locations = locationsAll,
    #                                fun = any)
    #  detectedByAll01 = rep(0, nPlumes(detectableSim))
    #  detectedByAll01[detectedByAll[[2]]] = 1
    #} else {
    detectedByAll = summaryPlumes(simulations = detectableSim,
                                  locations = locationsAll,
                                  values = detectable,
                                  fun = max)[[2]]
    #}
    nUndetectedAll = sum(1 - detectedByAll)
    if (nUndetectedAll > 0){
      paste0(nUndetectedAll, "plume(s) cannot be detected at any of the 'locationsAll' and is/are ignored.")
    }
  }


    # delete detected plumes, keep only locationsAll
    #if (dataType(detectableSim@values) == "LOG1S"){
    #  detectableFinal = subset(x = detectableSim,
    #                           plumes = which(detectedByFix01 == 0 & detectedByAll01 == 1),
    #                           locations = locationsAll)
    #} else {
      keepPlumes = which(detectedByFix == 0 & detectedByAll == 1)
      detectableFinal = subset(x = detectableSim,
                               plumes = keepPlumes,
                               locations = locationsAll)
    #}
  #} else {
  #  detectableFinal = subset(x = detectableSim,
  #                           locations = locationsAll)
  #}

  # can data be processed in memory
  if (!canProcessInMemory(detectableFinal@values, n = 2)){ # which n is realistic (how many copies needed in theory?)
    stop("Global search cannot be applied as values of simulations are too big to be loaded into memory.")
  }

  # extract values into matrix
  #if (dataType(detectableSim@values) == "LOG1S"){
  #  detectable_ = rep(0, prod(dim(detectableFinal@values)))
  #  detectable_[getValues(detectableFinal@values)] = 1
  #} else {
    detectable_ = getValues(detectableFinal@values)
  #}
  detectable_ = as.logical(detectable_)
  #print(table(detectable_))
  detectable = matrix(detectable_, byrow = TRUE, nrow = nrow(detectableFinal@values))

  # --------------- optimise ----------------- #
  # if detectable has 0 or 1 columns/rows

  if (dim(detectable)[1] == 1){
    if (dim(detectable)[2] == 1){
      result = list()
      result[["SD"]] = matrix(locationsAll[1])
      result[["cost"]] = 1
      result[["report"]] = list()
      result[["plumes"]] = keepPlumes
      result[["detectable"]] = detectable
      warning(paste0())
    } else {
      # setting dimnames only works, if first the dimension with extent of more than 1 is set
      dimnames(detectable)[[2]] = 1:dim(detectable)[2]
    }
  }
  dimnames(detectable)[[1]] = 1:dim(detectable)[1]
  dimnames(detectable)[[2]] = 1:dim(detectable)[2]

  # find optimal detection and one design to fulfil it
  #nameSave = NA
  if (!is.na(nameSave)){
    nameSaveF = paste0(nameSave, "_first")
  } else {
    nameSaveF = NA
  }
  firstOptimalDesign = completeSearch(Detectable = detectable,                              # matrix of locations(rows) and plumes(columns) of 0 - plume here not detectable and 1; rows not necessarily
                                      n = aimNumberNF,                                        # number of sensors
                                      lowerLimit = 1,                                       # only search for sensors that detect at least so many plumes
                                      increaseLimit = TRUE,                                 # if TRUE, the lower limit is increased to the current best cost + 1 every time the current best cost >= current lowerLimit; if FALSE lowerLimit is constant
                                      iterations = maxIterations,                           # maximal number of iterations, if no value given: number of all possibilities
                                      nameSave = nameSaveF                                   # path to save results, if it is a directory, it has to end with /
  )

  testLessSensors = FALSE
  if (firstOptimalDesign$lowerLimit > dim(detectable)[2]){ # "first" can detect all plumes -> check if even less sensors can do that
    if (findAllOptima){# to find all optima, their length must be known beforehand
      findSensorNumber = TRUE
    }
    if (findSensorNumber){
      testLessSensors = TRUE
      firstOptimalDesignsLessSensors = list()
      firstOptimalDesignsLessSensors[[1]] = firstOptimalDesign
      i = 1
      allDetected = TRUE
      while (i < aimNumberNF & allDetected){
        currentAimNumber = aimNumberNF - i
        if (!is.na(nameSave)){
          nameSaveFi = paste0(nameSave, "_first_", i)
        } else {
          nameSaveFi = NA
        }
        firstOptimalDesignsLessSensors[[i + 1]] =
          completeSearch(Detectable = detectable,                              # matrix of locations(rows) and plumes(columns) of 0 - plume here not detectable and 1; rows not necessarily
                         n = currentAimNumber,                                  # number of sensors
                         lowerLimit = 1,                                       # only search for sensors that detect at least so many plumes
                         increaseLimit = TRUE,                                 # if TRUE, the lower limit is increased to the current best cost + 1 every time the current best cost >= current lowerLimit; if FALSE lowerLimit is constant
                         iterations = maxIterations,                           # maximal number of iterations, if no value given: number of all possibilities
                         nameSave = nameSaveFi                                  # path to save results, if it is a directory, it has to end with /
          )

        allDetected = firstOptimalDesignsLessSensors[[i + 1]]$lowerLimit > dim(detectable)[2]
        i = i + 1
      }
      if (aimNumberNF == 1){ # then, no further 'firstOptimalDesignsLessSensors' are generated
        currentAimNumber = 1
      }
      # go back to last currentAimNumber when all plumes are detected
      detectionsByLessSensors = sapply(X = firstOptimalDesignsLessSensors, FUN = "[[", "lowerLimit")
      minAimNumberAllDetected = min((aimNumberNF:currentAimNumber)[detectionsByLessSensors > dim(detectable)[2]])

     warning(paste0(minAimNumberAllDetected, " sensors are sufficient to detect all plumes. Therefore optimal designs of this size are returned."))
     aimNumberNF = minAimNumberAllDetected
     firstOptimalDesign = firstOptimalDesignsLessSensors[[max(which(detectionsByLessSensors > dim(detectable)[2]))]]
    }
  }
  if (findAllOptima){
  # find all optimal designs
    if (!is.na(nameSave)){
      nameSaveA = paste0(nameSave, "_all")
    } else {
      nameSaveA = NA
    }
    allOptimalDesign =  completeSearch(Detectable = detectable,                              # matrix of locations(rows) and plumes(columns) of 0 - plume here not detectable and 1; rows not necessarily
                                       n = aimNumberNF,                                        # number of sensors
                                       lowerLimit = firstOptimalDesign[["lowerLimit"]] - 1,                                       # only search for sensors that detect at least so many plumes
                                       increaseLimit = FALSE,                                # if TRUE, the lower limit is increased to the current best cost + 1 every time the current best cost >= current lowerLimit; if FALSE lowerLimit is constant
                                       iterations = maxIterations,                           # maximal number of iterations, if no value given: number of all possibilities
                                       nameSave = nameSaveA                                  # path to save results, if it is a directory, it has to end with /
    )
  }

  # ------------- reconstruct optimal designs ignored by merging ------------ #
  if (!findAllOptima){# find one optimum (equal SDs from merging are not reconstructed)
    SD = as.integer(firstOptimalDesign$finalSD[nrow(firstOptimalDesign$finalSD),])
    if (is.element("0", firstOptimalDesign$finalSD)){
      SD = SD[SD != 0]
    }
  } else {
    # find all optima (reconstruct all equal SDs)
    # reconstruct sensors that detect the same remaining plumes
    SDsensors = list()
    for(k in 1:nrow(allOptimalDesign$finalSD)){
      SDsensors[[k]] = list()
      for(l in 1:aimNumberNF){
        if(l == 1){
          alreadyDetected = rep(FALSE, dim(detectable)[2])
        }else{
          alreadyDetected = apply(X = detectable[as.integer(allOptimalDesign$finalSD[k,1:(l-1)]),,drop = FALSE], # changed: insterted "as.integer"
                                  MARGIN = 2, FUN = any)
        }
        currentDetectable = detectable[,!alreadyDetected,drop = FALSE]
        hereDetected = currentDetectable[as.integer(allOptimalDesign$finalSD[k,l]),,drop = FALSE] # changed: insterted "as.integer"
        whichSame = which(
          apply(
            X = matrix(
              apply(X = currentDetectable, MARGIN = 1, FUN = "==", hereDetected),
              nrow = length(hereDetected)),
            MARGIN = 2, FUN = all)
        )
        SDsensors[[k]][[l]] = whichSame #as.integer(dimnames(currentDetectable)[[1]][whichSame])
      }
    }
    if (!is.na(nameSave)){
      save(SDsensors, file = paste0(nameSave, "_SDsensors.Rdata"))
      if (verbatim){
        print(paste0("All locations that are equivalent sensors have been determined and saved at ", nameSave, "_SDsensors.Rdata"))
      }
    }

    # generate all relevant sensor combinations
    ## number of combination
    nCand = matrix(nrow = length(SDsensors), ncol = aimNumberNF)
    for(k in seq(along = SDsensors)){
      nCand[k,] = sapply(SDsensors[[k]], length)
    }

    SDmatrix = matrix(nrow = sum(apply(nCand, 1, prod)), ncol = aimNumberNF)
    m = 0
    for (k in seq(along = nCand[,1])){
      # 1st sensors
      SDmatrix[m + 1:prod(nCand[k,]), 1] = rep(SDsensors[[k]][[1]], each = prod(nCand[k,-1]))
      if (aimNumberNF > 2){
        # 2nd - (aimNumberNF -1)st sensors
        for(l in 2:(aimNumberNF - 1)){
          SDmatrix[m + 1:prod(nCand[k,]), l] =
            rep(
              rep(SDsensors[[k]][[l]], each = prod(nCand[k,(l + 1):aimNumberNF])),
              times = prod(nCand[1:(l-1)]))
        }
      }
      # aimNumberNF-st sensor
      SDmatrix[m + 1:prod(nCand[k,]), aimNumberNF] =
        rep(SDsensors[[k]][[aimNumberNF]], prod(nCand[k,-aimNumberNF]))
      m = m + prod(nCand[k,])
    }
    if (!is.matrix(SDmatrix)){
      SDmatrix = matrix(SDmatrix, nrow = 1)
    }
    # make them unique (when including equal locations, some of them may have been used earlier and should not be in; they need deletion afterwards)
    sorted = t(apply(SDmatrix, 1, sort))
    duplicated = duplicated(sorted)
    SDmatrix = SDmatrix[!duplicated,, drop = FALSE]
    if (!is.na(nameSave)){
      save(SDmatrix, file = paste0(nameSave, "_SDmatrix.Rdata"))
      if (verbatim){
        print(paste0("All equivalent sampling designs have been determined and saved at ", nameSave, "_SDmatrix.Rdata"))
      }
    }

    # test if all optimal (they all should be it anyway)
    nDetections = integer(nrow(SDmatrix))
    for(k in seq(along = SDmatrix[,1])){
      nDetections[k] =
        sum(apply(detectable[SDmatrix[k,],, drop = FALSE], 2, any))
    }
    isOptimal = nDetections == allOptimalDesign[["lowerLimit"]]
    if (!all(isOptimal)){
      warning(paste0("Some of the determined sampling designs are not optimal (this should not happen). They are deleted."))
      SDmatrixOpt = SDmatrix[isOptimal, , drop = FALSE]
    } else {
      SDmatrixOpt = SDmatrix
    }
   if (!is.na(nameSave)){
     save(SDmatrix, SDmatrixOpt, file = paste0(nameSave,"_SDmatrix.Rdata"))
   }
  }
  # until here locations in SD/SDmatrix refer to detectable as used as input to completeSearch,
  # it ignores previous deleting of locations (fix, all)
  if (!findAllOptima){
    SD = matrix(locationsAll[SD], nrow = 1)
  } else {
    SD = matrix(locationsAll[SDmatrixOpt], ncol = aimNumberNF)
  }
  SDf = matrix(nrow = nrow(SD), ncol = ncol(SD) + length(locationsFix))
    for (i in seq(along = SD[,1])){
      SDf[i,1:ncol(SD)] = SD[i,]
      if (length(locationsFix) > 0){
        SDf[i,ncol(SD) + 1:length(locationsFix)] = locationsFix
      }
      SDf[i,] = sort(SDf[i,])
    }

  # cost: fraction of missed plumes
  simulationsAtSensors = subset(detectableSim, locations = c(SD[1,], locationsFix))
  cost = (nPlumes(simulations) - sum(apply(matrix(getValues(simulationsAtSensors@values),
                                                  byrow = TRUE, ncol = nPlumes(simulations)),
                                           2, any)))/nPlumes(simulations)
  evaluation = data.frame(cost = cost, number = ncol(SD))
  # ---------------- output ------------------------------------ #
  result = list()
  result[["SD"]] = list(SDf)
  result[["evaluation"]] = evaluation
  result[["report"]] = list()
  result[["report"]][["detectable"]] = detectable
  result[["report"]][["first"]] = firstOptimalDesign
  if (testLessSensors){
    result[["report"]][["SDsLessSensors"]] = firstOptimalDesignsLessSensors
  }
  if (findAllOptima){
    result[["report"]][["all"]] = allOptimalDesign
  }
#   if (findAllOptima){
#     result[["SD"]] = SD#allOptimalDesign[["finalSD"]]
#     result[["cost"]] = cost # allOptimalDesign[["lowerLimit"]]
#     result[["report"]] = list()
#     result[["report"]][["first"]] = firstOptimalDesign
#
#   } else {
#     result[["SD"]] = SD #firstOptimalDesign[["finalSD"]][length(firstOptimalDesign[["finalSD"]]),] # last design is the optimal one
#     result[["cost"]] = #firstOptimalDesign[["lowerLimit"]] - 1
#     result[["report"]] = list()
#     result[["report"]][["first"]] = firstOptimalDesign
#   }
#   result[["plumes"]] = keepPlumes
#   result[["detectable"]] = detectable
#   if (testLessSensors) result[["report"]][["SDsLessSensors"]] = firstOptimalDesignsLessSensors
  return(result)
}
