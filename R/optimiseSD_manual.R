############################################################
# optimiseSD_manual                                        #
############################################################
# generate map

optimiseSD_manual = function(simulations,
                             costFun,
                             locationsAll = 1:nLocations(simulations),
                             locationsFix = integer(0),
                             locationsInitial = integer(0),
                             aimCost = NA, aimNumber = NA,
                             nameSave = NA, plot = TRUE, verbatim = FALSE,
                             costMap = NA,
                             maxIterations = 10,
                             valuesPlot = integer(0), # locations$data to plot in addition to the computed cost map (by default none)
                             colors = grey.colors
                             ){

  # ------------------------- prepare  ------------------------------ #
#  if (is.na(costMap)){
#    stop("'costMap' required.")
#  }
  # - - - - clean location indices - - - - - - -
  locations = cleanIndices(
    locationsTotal = 1:nLocations(simulations),#length(simulations@locations),
    locationsAll = locationsAll,
    locationsFix = locationsFix,
    locationsInitial = locationsInitial
  )
  locationsAll = locations[["locationsAll"]]
  locationsFix = locations[["locationsFix"]]
  locationsCurrent = locations[["locationsInitial"]]

  # - - - - - initialise indices of change - - - - - -
  locationsAdd = integer(0)
  locationsDel = integer(0)
  addORdel = ""

  # - - - - - define and initialise report - - - - -
  report = list()
  report[["cost"]] = numeric(maxIterations + 1)
  report[["identify"]] = list()
  report[["locationsCurrent"]] = list()

  report[["identify"]][[1]] = list()
  report[["locationsCurrent"]][[1]] = locationsCurrent

  # - - - -  set plotting parameters - - - - - #
  # points for plotting
  if (class(simulations@locations) == "SpatialIndexDataFrame"){
    simPoints = data.frame(y = 1:nLocations(simulations), x = 1)
    coordinates(simPoints) = ~ x + y
  } else {
    simPoints = SpatialPoints(
      coords = coordinates(simulations@locations),
      proj4string = CRS(proj4string(simulations@locations)))
  }


  # point settings: all, fix, current; added, deleted; invalid  (both: in last step)
  colPoints = c(rep(4, 3), 3, 2, 7)
  pchPoints = c(8, 20, 20, 20, 1, 20)
  cexPoints = c(1, 0.1, 1, 1, 1, 0.1)

  # trellis settings
  settings = list(regions = list(col = colors),
                  superpose.symbol = list(
                    col =  colPoints,
                    pch = pchPoints,
                    cex = cexPoints))
  trellis.par.set(settings)


  for (i in 1:maxIterations){
    #  - - - - - cost map - - - - - - - -
    # compute
    map = costMap(
      simulations = simulations,
      locations = c(locationsFix, locationsCurrent)
      #nameSave = nameSave
    )
    report[["cost"]][i] = map[["cost"]]

    # attach to spatial object for plotting
    simulations@locations@data$costMap = map[["costLocations"]]


    # - - - - - points - - -- - - - - -
    # current points for plotting
    plotPoints = list(list('sp.points', simPoints[locationsFix,], col = colPoints[1], pch = pchPoints[1], cex = cexPoints[1]),
                      list('sp.points', simPoints[locationsAll,], col = colPoints[2], pch = pchPoints[2], cex = cexPoints[2]),
                      list('sp.points', simPoints[locationsCurrent,], col = colPoints[3], pch = pchPoints[3], cex = cexPoints[3]))

    if (addORdel == "a"){
      plotPoints[[4]] =  list('sp.points', simPoints[locationsAdd,], col = colPoints[4], pch = pchPoints[4], cex = cexPoints[4])
    }
    if (addORdel == "d"){
      plotPoints[[4]] =  list('sp.points', simPoints[locationsDel,], col = colPoints[5], pch = pchPoints[5], cex = cexPoints[5])
    }

    # --------------------- plot ----------------------------------------

    # plot cost map with locations
    if (class(simulations@locations) == "SpatialPolygridDataFrame"){ # plotted SDF is not original ->
      # to find original from identify in plot, reconstruction must be possible
      mapPlot0 <- spplot(simulations@locations, zcol = c("costMap", valuesPlot),
                         sp.layout = plotPoints, col.regions = colors,
                         main = paste0("cost: ", map[[1]]),
                         returnSGDF = TRUE)
      mapPlot = mapPlot0[["spplot"]]
      coordIndex = mapPlot0[["grid"]]@data$index
    } else { # plotted SDF is original (?)
      mapPlot <- spplot(simulations@locations, zcol = c("costMap", valuesPlot),
                        sp.layout = plotPoints, col.regions = colors,
                        main = paste0("cost: ", map[[1]]))
      coordIndex = 1:nLocations(simulations) # does this work for all kind of SCD except SPolygridDF
    }

    # add point key at bottom + space for key
    mapPlot1 = update(mapPlot,
                      zcol = c("costMap", valuesPlot),
                      key = simpleKey(
                        c("fix, cannot be changed",
                          "potential, may be added",
                          "current, may be deleted",
                          "added in last iteration",
                          "deleted in last iteration"),
                        points = TRUE, columns = 1,
                        space = "bottom"))

    plot(mapPlot1)

    # interactive
    addORdel = ""
    while (!is.element(addORdel, c("a", "d", "s"))){
      addORdel = readline("Add ('a') or delete ('d') sensors, or stop ('s')?")
    }
    if (addORdel != "s"){
      trellis.focus()
      selected = panel.identify()
      if (class(simulations@locations) == "SpatialPolygridDataFrame"){
        selected1 = unique(coordIndex[selected])
        if (length(selected1) < length(selected)){
          print(paste0(length(selected1) - length(selected1), "of the selected locations are actually identical to other selected ones, they are ignored."))
        }
      } else {
        selected1 = selected
      }

      if (addORdel == "a"){
        # determine selected locations that can actually be added
        locationsAdd = intersect(selected1, locationsAll)

        if (length(locationsAdd) < length(selected1)){
          print(paste0(length(selected1) - length(locationsAdd), " selected locations are ignored as they are not part of 'locationsAll'."))
        }
        # update locations
        locationsAll = setdiff(locationsAll, locationsAdd)
        locationsCurrent = union(locationsCurrent, locationsAdd)
      }
      if (addORdel == "d"){
        # determine selected locations that can actually be added
        locationsDel = intersect(selected1, locationsCurrent)
        if (length(locationsDel) < length(selected1)){
          print(paste0(length(selected1) - length(locationsDel), " selected locations are ignored as they were not part of the current locations."))
        }

        # update locations
        locationsAll = union(locationsAll, locationsDel)
        locationsCurrent = setdiff(locationsCurrent, locationsDel)
      }
    } else { # addORdel = "s"
      break
    }
    report[["identify"]][[i + 1]] = selected
    report[["locationsCurrent"]][[i + 1]] = locationsCurrent
  }
  report[["cost"]] = report[["cost"]][1:i]


  SDraw =list()
  for (i in 1:i){
    SDraw[[i]] = sort(c(report[["locationsCurrent"]][[i]], locationsFix))
  }
  nSDsAll = sapply(X = SDraw, FUN = length)
  nSDs = sort(unique(nSDsAll))
  SD = list()
  eval = data.frame(cost = numeric(length(nSDs)), number = integer(length(nSDs)))
  for (i in seq(along = nSDs)){
    whichNi = which(nSDsAll == nSDs[i])
    cost_i = numeric(length(whichNi))
    for (j in seq(along = whichNi)){
      cost_i[j] = costMap(
        simulations = simulations,
        locations = SDraw[[whichNi[j]]]
      )[["cost"]]
    }
    # extract the best
    best_i = which(cost_i == min(cost_i))
    SD[[i]] = matrix(nrow = length(best_i), ncol = nSDs[i])
    for (j in seq(along = best_i)){
      SD[[i]][j,] = SDraw[[whichNi[best_i[j]]]]
    }
    SD[[i]] = unique(SD[[i]])
    eval$cost[i] = min(cost_i)
    eval$number[i] = nSDs[i]
  }
  out = list()
  out[["SD"]] = SD
  out[["evaluation"]] = eval
  out[["report"]] = report
  return(out)
}
