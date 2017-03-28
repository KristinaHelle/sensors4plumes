# to be decided:
# a) ways to define aim as nr of sensors: aimN, length of locationsInitial?
# b) how to treat multiple/impossible initial or possible locations: delete with warning?
# c) stopping rules depend on algorithm

# a) ways to define initial
#   greedy: if not explicitly defined, use 0 sensors and add with algorithm
#   ssa: either define explicitly OR give number and use simple random (for spatially balanced random etc. define explicitly beforehand with other fun)
#   gen: rbga.bin can take several initial sampling designs (each in row of matrix) or do it randomly
#   complete: cannot take it into account, may use it to define number
## initial: can be vector or list; if list given and not gen used, only first entry is taken into account
# b) delete with warning
# c) define all parameters specific for algorithm beforehand (via replaceDefault)

# change: only use concrete default values, not names of values to be impoarted from wrapping function; for import define there what to use


cleanIndices = function(
  locationsTotal,
  locationsAll = locationsTotal,
  locationsFix = integer(0),
  locationsInitial = integer(0)
){
  # kind of initial locations
  multipleInitial = FALSE
  if (is(locationsInitial, "matrix")){
    multipleInitial = TRUE
    nInitial = ncol(locationsInitial)
  } else {
    nInitial = length(locationsInitial)
  }

  # all locations must exist
  if (!all(is.element(c(locationsAll, locationsFix, locationsInitial), locationsTotal))) {
    warning("Some of the 'locations' do not exist, they are ignored.")

    locationsAll = intersect(locationsAll, locationsTotal)
    locationsFix = intersect(locationsFix, locationsTotal)

    if (multipleInitial){
      for (i in 1:dim(locationsInitial)[1]){
        locationsInitial[i,!is.element(locationsInitial[i,], locationsTotal)] = NA
      }
    } else {
      locationsInitial = intersect(locationsInitial, locationsTotal)
    }
  }

  # locationsAll must contain all locationsInitial
  if (!all(is.element(locationsInitial, locationsAll))){
    warning("Some of the 'locationsInitial' are not part of 'locationsAll'. The union of both sets is used as possible sensor locations (instead of 'locationsAll').")
    locationsAll = union(locationsAll, na.omit(as.integer(locationsInitial)))
  }

  # locationsFix must not be part of either locationsInitial or locationsAll
  if (any(is.element(locationsFix, c(locationsAll)))){
    warning("'locationsFix' overlaps with 'locationsAll' or 'locationsInitial'; overlapping locations are regarded as fix.")
    locationsAll = setdiff(locationsAll, locationsFix)
    if (multipleInitial){
      for (i in 1:dim(locationsInitial)[1]){
        locationsInitial[i, is.element(locationsInitial[i,], locationsFix)] = NA
      }
    } else {
      locationsInitial = setdiff(locationsInitial, locationsFix)
    }
  }

  # delete multiples in initial locations (for a matrix this cannot be done by 'unique')
  if (multipleInitial){
    for (i in 1:dim(locationsInitial)[1]){
     locationsInitial[i, duplicated(locationsInitial[i,])] = NA
    }
  } else {
    locationsInitial = unique(locationsInitial)
  }
  # delete invalid initial locations
  if (multipleInitial){
    nNAinitial = apply(FUN = sum, X = is.na(locationsInitial), MARGIN = 1)
    if (min(nNAinitial) == nInitial){# all NA
      warning("No valid 'locationsInitial'.")
      locationsInitial = integer(0)
    } else {
      if (length(unique(nNAinitial)) == 1){# all rows have same number of valid entries
        new_locationsInitial = matrix(nrow = nrow(locationsInitial), ncol = nInitial - nNAinitial)
        for (i in 1:dim(locationsInitial)[1]){
          new_locationsInitial[i,] = na.omit(locationsInitial[i,])
        }
      } else {# keep rows with most valid entries
        longestInitial = which(nNAinitial == min(nNAinitial))
        new_locationsInitial = matrix(nrow = length(longestInitial), ncol = nInitial - min(nNAinitial))
        for (i in seq(along = longestInitial)){
          new_locationsInitial[i,] = na.omit(locationsInitial[longestInitial[i],])
        }
        warning(paste0("Some 'locationsInitial' had less valid entries, they are ignored: ",
                       length(longestInitial), " row(s) of 'locationsInitial' left."))
      }
      locationsInitial = new_locationsInitial
    }
  }
  if (multipleInitial & length(locationsInitial) > 0){
    if (ncol(locationsInitial) < nInitial){
      warning(paste0("Each row of 'locationsInitial' had invalid entries: ", ncol(locationsInitial), " column(s) left."))
    }
  } else {
    if (length(locationsInitial) < nInitial){
      warning(paste0("Some entries of 'locationsInitial' were invalid: ", length(locationsInitial), " entry(es) left."))
    }
  }

  locations = list()
  locations[["locationsAll"]] = unique(setdiff(union(locationsAll, locationsInitial), locationsFix))
  locations[["locationsFix"]] = unique(locationsFix)
  locations[["locationsInitial"]] = locationsInitial
  #locations[["locationsInitial"]] = unique(setdiff(locationsInitial, locationsFix))

  return(locations)
}

optimiseSD = function(
  simulations,                                      # needed? yes, defines possible locations, values to be used in cost function
  costFun,
  locationsAll = 1:nLocations(simulations),
  locationsFix = integer(0),
  locationsInitial = integer(0),                    # can be a list (for genetic)
  aimCost = NA,
  aimNumber = NA,
  optimisationFun,
  nameSave = NA, #"optimiseSD", # to be used as part of the filenames; kind of files and filenames may depend on algorithm
  plot = FALSE,
  verbatim = FALSE,
  ...
  ){
  costFun = replaceDefault(costFun, newDefaults = list(simulations = simulations))[[1]]
  # ------------------------- prepare location indices ------------------------------ #
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

  # -------------------- optimisation --------------------------------------- #
  #costInitial = costFun(simulations = simulations,
  #                      locations = locationsInitial)
                        #plot = FALSE)
  optimalSD = optimisationFun(simulations = simulations,
                                     costFun = costFun,
                                     locationsAll = locationsAll,
                                     locationsFix = locationsFix,
                                     locationsInitial = locationsInitial,
                                     aimCost = aimCost,
                                     aimNumber = aimNumber,
                                     nameSave = nameSave,
                                     verbatim = verbatim)

  # determine which SD is wanted by aimCost or aimNumber
  aim = list()
  if (!is.na(aimNumber)){
    aim[["SDaimNumber"]] = which(optimalSD$evaluation$number == aimNumber)
  }
  if (!is.na(aimCost)){
    whichBelow = which(optimalSD$evaluation$cost <= aimCost)
    if (length(whichBelow) > 1){
      aim[["SDaimCost"]] = min(whichBelow)
    } else {
      aim[["SDaimCost"]] = whichBelow
    }
  }
  optSD = list(
    SD = optimalSD[["SD"]],
    evaluation = optimalSD[["evaluation"]],
    aimSD = aim,
    report = optimalSD[["report"]]
  )

  return(optSD)
}
  # = do.call(what = optimisationFun, args = as.list(formals(optimisationFun)))


#   optimalSD = optimisationFun(
#     simulations = simulations,                     # simulations object to optimise sensors for; forwarded to 'costFun', 'locations...' must be subset of locations else they are ignored
#     costFun = costFun,                         # function to compute cost
#     locationsAll = locationsAll,                    # all possible sensor locations (as indices of simulations@locations)
#     locationsFix = locationsFix,                    # these sensors are always included to determine cost, they cannot be deleted
#     locationsInitial = locationsInitial,                # current sensors that can be deleted
#     aimCost = aimCost,
#     aimNumber = aimNumber,
#     nameSave = nameSave                       # without suffix .Rdata
#     #  plot = FALSE,
#   )


  #costI = costInitial[[1]]
  #whichAim = c(!missing(aimNumber), !missing(aimCost)) # CHANGE
#   switch(algorithm,
#     "SSA" = {
#       optimalSD = optimiseSD_ssa(
#         simulations = simulations,
#         costFun = costFun,
#         locationsAll = locationsAll,
#         locationsFix = locationsFix,
#         locationsInitial = locationsInitial,
#         aim = aimCost, # CHANGE
#         plot = plot,
#         verbatim = verbatim,
#         pathSave = pathSave,
#         nameSave = nameSave
#         )
#     },
#     "greedy" = {
#       # determine if optimisation aims at a maximal number of sensors or at a cost limit
#
#       # run greedy optimisation
#       optimalSD = optimiseSD_greedy (
#         simulations = simulations,
#         costFun = costFun,
#         locationsAll = locationsAll,
#         locationsFix = locationsFix,
#         locationsInitial = locationsInitial,
#         aimCost = aimCost,
#         aimNumber = aimNumber,
#         nameSave = nameSave #paste0(pathSave, "/", nameSave)
#         #plot = plot
#       )
#
#       # check secondary aim if required
#       if (identical(whichAim, c(TRUE, TRUE))){
#         if(optimalSD[["evalSDs"]]$cost[optimalSD[["finalSDwhich"]][1]] <= aimCost){# CHANGE
#          message("The desired cost was achieved as well.")
#         }else{
#           message(paste("The desired cost was not achieved, cost of optimal SD is ", optimalSD[["evalSDs"]]$cost[optimalSD[["finalSDwhich"]][1]], sep = ""))
#         }
#       }
#     },
#     "genetic" = {
#       require(genalg)
#       cost0sensors = costFun(simulations, locations = integer(0))
#       # chromosome must contain only the locations that can be changed, not fix, not forbidden
#       # suggestions, from locationsInitial
#       if (!missing(locationsInitial)){
#         suggestionsTotal = rep(0, nLocations(simulations))
#         suggestionsTotal[c(locationsInitial, locationsFix)] = 1
#         suggestions = matrix(suggestionsTotal[locationsAll], nrow = 1)
#       }
#       # prepare evalFunc
#       # what to do if no aimNumber / aimCost given?
#       evalFunc = function(chromosome){
#         out = evalFun(
#           chromosome = chromosome,
#           simulations = simulations,
#           costFun = costFun,
#           locationsAll = locationsAll,
#           locationsFix = locationsFix,
#           aimNumber = aimNumber,#CHANGE
#           aimCost = aimCost,# CHANGE
#           plot = FALSE)[[1]]
#       }
#       monitorFunc = function(rbgaResults){
#         out = plotSD(
#           simulations = simulations,
#           costFun = costFun,
#           type = "genetic",
#           SD = locationsOld,
#           plot = plot,
#           param = list(
#             locationsFix = locationsFix,
#             locationsInitial = locationsInitial,
#             locationsAll = locationsAll,
#             rbgaResults = rbgaResults
#           )
#         )
#       }
#         if (missing(aimNumber)){
#           aimNumber = length(unique(c(locationsInitial, locationsFix)))
#         }
#         nAll = length(locationsAll)
#         nFix = length(locationsFix)
#         rbgaResults = rbga.bin(size = nAll,
#                                 suggestions = suggestions,
#                                 evalFunc = evalFunc,
#                                 zeroToOneRatio = nAll/(aimNumber - nFix),
#                                 monitorFunc = monitorFunc,
#                                 verbose = verbatim)
#
#       optimalSD = SDrbga(simulations = simulations,
#                 locationsFix = locationsFix,
#                 locationsAll = locationsAll,
#                 rbgaResults = rbgaResults,
#                 costFun = costFun)
#       optimalSD[["report"]] = rbgaResults
#       # aimNumber: penalty, if too much
#       # aimCost: stop if reached
#    },
#     "global" = {
#       # detectable = simulatios@values must be in memory
#       # and consist of 0 and  1 only
#
#       completeSearch(Detectable = detectable,                              # matrix of locations(rows) and plumes(columns) of 0 - plume here not detectable and 1; rows not necessarily
#                               n = length(locationsInitial))
#
#     }
#   )
