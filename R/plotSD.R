
# it would be nice to have logarithmic scale with according key
plotSD = function(
  simulations,
  SD,
  locationsFix = integer(0),
  locationsInitial = integer(0),
  locationsAll = integer(0),
  costMap,
  zcol = 1, # background (from locations@data); not used if cost map given
  allIn1Plot = 0, # 0: not all plots in one, >0: all in one with this SDs map as background
  pch = c(1, 20, 4),
  col = c("white", "white", "white", "white"),
  pointsKey = TRUE,
  mainCost = TRUE,
  pch.SDs, # if each SD is to be plotted with different style
  col.SDs,
  cex.SDs,
  ...
#  layout # to be forwarded to sp.layout (in addition to points from SD etc.)
){
  # zcol must be character
  if (!is.character(zcol)){
    zcol = names(simulations@locations@data)[zcol]
  }
  
  # points at simulations locations
  if (class(simulations@locations) == "SpatialIndexDataFrame"){
    simPoints = data.frame(y = 1:nLocations(simulations), x = 1)
    coordinates(simPoints) = ~ x + y
  } else {
    simPoints = SpatialPoints(
      coords = coordinates(simulations@locations), 
      proj4string = CRS(proj4string(simulations@locations)))    
  }
  
  # SD: turn into list
  if (!is.list(SD)){
    if (is.matrix(SD)){
      SDs = list()
      for (i in 1:nrow(SD)){
        SDs[[i]] = SD[i,]
      }
      SD = SDs
    } else {# single SD
      SD = list(SD)      
    }
  }
  
  # plot SDs as circles of different size or as customised
  plotSDs = list()
  cexFactor = 1/length(SD)
  if (missing(pch.SDs)){
    pch.SDs = rep(pch[1], length(SD))
  }
  if (missing(col.SDs)){
    col.SDs = rep(col[1], length(SD))
  }
  if (missing(cex.SDs)){
    cex.SDs = 1/length(SD) * 1:length(SD)
  }
  for (i in seq(along = SD)){
    plotSDs[[i]] = list("sp.points", simPoints[SD[[i]],], 
                        col = col.SDs[i], 
                        pch = pch.SDs[i], 
                        cex = cex.SDs[[i]])
  }   
  plotFix = list("sp.points", simPoints[locationsFix,], 
                    col = col[2], pch = pch[2])
  plotInitial = list("sp.points", simPoints[locationsInitial,], 
                    col = col[3], pch = pch[3])
  plotAll = list("sp.points", simPoints[locationsAll,], 
                 col = col[4], pch = ".")
  cex = c(1,1,1,0.2)
  pch = c(pch, 15)
  col[col == "white"] = 1
  settings = list(superpose.symbol = list(
                    col =  col, 
                    pch = pch,
                    cex = cex))  
  trellis.par.set(settings)
  
  if (missing(costMap)){
    # plot with given data as background
    plots = spplot(simulations@locations, zcol = zcol,
                      sp.layout = list(plotSDs, plotAll, plotFix, plotInitial),
                      ...)
    if (pointsKey){
      plots = update(plots, 
                     key = simpleKey(
                       c("sampling design(s)",
                         "fix sensors",
                         "initial sensors",
                         "potential sensor locations"),
                       points = TRUE, columns = 1, 
                       space = "bottom"))      
    }
  } else {
    if (allIn1Plot == 0){# not all plots in one
      plots = list()
      for (i in seq(along = SD)){
        # compute cost map(s)
        computeCostMap = costMap(simulations = simulations,
                                 locations = c(SD[[i]], locationsFix))
        if (!is.element("costLocations", names(computeCostMap))){
          stop("Cost map could not be derived, output of 'costMap' does not contain 'costLocations'.")
        }
        # generate plot(s)
        simulations@locations@data[,paste0("costMap_",i)] = computeCostMap[["costLocations"]]
        plots0 = spplot(simulations@locations, 
                            zcol = paste0("costMap_",i),
                            sp.layout = list(plotSDs[[i]], plotAll, plotFix, plotInitial),
                            ...)   
        if (is.element("cost", names(computeCostMap)) & mainCost){
          plots0 = update(plots0, main = signif(computeCostMap[["cost"]], 5))
        }
        if (pointsKey){
          plots[[i]] = update(plots0, 
                              key = simpleKey(
                                c("sampling design(s)",
                                  "fix sensors",
                                  "initial sensors",
                                  "potential sensor locations"),
                                points = TRUE, columns = 1, 
                                space = "bottom"))          
        } else {
          plots[[i]] = plots0
        }
      }      
    }
    if (allIn1Plot > 0){
      computeCostMap = costMap(simulations = simulations,
                               locations = c(SD[[allIn1Plot]], locationsFix))
      simulations@locations@data[,paste0("costMap_", allIn1Plot)] = 
        computeCostMap[["costLocations"]]
      
      plots = spplot(simulations@locations, 
                          zcol = paste0("costMap_",allIn1Plot),
                          sp.layout = list(plotSDs, plotAll, plotFix, plotInitial),
                          ...)
      if (is.element("cost", names(computeCostMap)) & mainCost){
        plots = update(plots, main = signif(computeCostMap[["cost"]], 5))
      }
      if (pointsKey){
        plots = update(plots, 
                       key = simpleKey(
                         c("sampling design(s)",
                           "fix sensors",
                           "initial sensors",
                           "potential sensor locations"),
                         points = TRUE, columns = 1, 
                         space = "bottom")) 
      }
    }
    return(plots)
  }
}







# plotSD = function(
#   simulations,
#   costFun,
#   type = NA, # ssa, greedy, genetic, global
#   SD,
#   plot = TRUE,
#   style = list(),
#   param = list(),
#   main = ""
#   ){
#   # 
#   #require(sp)
#   spatial = extractSpatialDataFrame(obj = simulations)
#   coord = coordinates(simulations@locations)
#   
#   # compute cost (and map)
#   if (is.na(type) | is.element(type, c("ssa", "greedy"))){
#     costCurrent = costFun(simulations = simulations, locations = SD, plot = plot)
#   }
#   if (!is.na(type)){
#     if (type == "genetic"){
#       optimalSD = SDrbga(simulations = simulations, 
#                          locationsFix = param[["locationsFix"]],
#                          locationsAll = param[["locationsAll"]],
#                          rbgaResults = param[["rbgaResults"]], 
#                          costFun = costFun)
#       print(optimalSD[["SDs"]])
#       SD1 = sample(optimalSD[["SDs"]], 1)[[1]]
#       costCurrent = costFun(simulations = simulations, locations = SD1, plot = plot)
#     }
#   }
#   
#   # transfer map values to spatial object
#   if (!is.null(costCurrent[["map"]])) {
#     spatial@data[,1] = costCurrent[["map"]]# rename to "map"
#   }
#   
#   # prepare other objects to be plotted 
#   if (is.element("locationsFix", names(param))){
#     style[["Fix"]] = list("sp.points", coord[
#       param[["locationsFix"]], , drop = FALSE], pch = 2, col = "grey")    
#   }
#   if (is.element("locationsInitial", names(param))){
#   style[["Init"]] = list("sp.points", coord[
#     param[["locationsInitial"]], , drop = FALSE],  pch = 1, col = "grey")
#   }
#   
#   if (is.na(type)){
#     style[["SD"]] = list("sp.points", coord[
#       SD,, drop = FALSE], pch = 20, col = "green")
#     main = paste("cost = ", signif(costCurrent[[1]], digits = 5), sep = "")
#   } else {
#     if (type == "ssa"){
#       style[["SD"]] = list("sp.points", coord[
#         SD,, drop = FALSE], pch = 20, col = "green")
#       style[["Potential"]] = list("sp.points", coord[
#         param[["locationsPotential"]], , drop = FALSE], pch = ".", col = "red")
#       style[["Change"]] = list("sp.points", coord[
#         param[["locationsChange"]], , drop = FALSE], pch = 4, col = "red")
#       
#       main = paste("cost = ", signif(costCurrent[[1]], digits = 5), "\n", 
#                    "best cost = ", signif(param[["costBest"]], digits = 5), "\n", 
#                    "Iterations = ", param[["iterations"]])
#     }  
#     
#     if (type == "genetic"){
#       style[["SD"]] = list("sp.points", coord[
#         SD1,, drop = FALSE], pch = 19, col = "green")
#       style[["SDs"]] = list("sp.points", coord[
#         unlist(optimalSD[["SDs"]]),, drop = FALSE], pch = 20, col = "green")
#       main = paste("cost = ", signif(costCurrent[[1]], digits = 5), "\n",
#                    "Iterations = ", param[["rbgaResults"]][["iter"]])
#     }  
#   }
#   spatial@data$total = apply(FUN = mean, X = as.matrix(spatial@data), MARGIN = 1, na.rm = TRUE)
# #  return(spatial)
#   # plot
#    print(spplot(spatial, zcol = "total",
#                 col.regions = bpy.colors(100),
# #                sp.layout = style,
#                 main = main))
#    
# }
# 
