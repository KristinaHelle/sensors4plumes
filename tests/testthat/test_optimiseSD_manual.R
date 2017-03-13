########################################################################
# test optimiseSD_manual                                               #
########################################################################

# test with all kinds of SpatialDataFrame in simulations@locations
# - - - - - - - - - SPolygridDF - - - - - - -
data(radioactivePlumes_local)
locAll1 = sample.int(nLocations(radioactivePlumes_local), 5000)
locFix1 = sample.int(nLocations(radioactivePlumes_local), 10)
locInit1 = sample.int(nLocations(radioactivePlumes_local), 10)
radioactivePlumes_local@locations@data$p1 = getValues(
  subset.Simulations(radioactivePlumes_local, plumes = 1, values = 1)@values)
meanFun = function(x){mean(x, na.rm = TRUE)}

map1 = spplot(radioactivePlumes_local@locations, zcol = "p1", returnSGDF = TRUE)
radioactivePlumes_local@locations@data$cost1 = spatialSpread(simulations = radioactivePlumes_local, 
                                                             locations = c(locFix1, locInit1), 
                                                             weightByArea = TRUE, 
                                                             fun = minimalDistance, 
                                                             fun_R = meanFun)[["costLocations"]]
map2 = spplot(radioactivePlumes_local@locations, zcol = "cost1", returnSGDF = TRUE)


# smaller
radioactivePlumes_loc = subset(radioactivePlumes_local, locations = 1:2001, overwrite = TRUE)
map3 = spplot(radioactivePlumes_loc@locations, zcol = "cost1", returnSGDF = TRUE)

spatialSpreadMinDist = replaceDefault(
  spatialSpread,
  newDefaults = list(
    weightByArea = TRUE, 
    fun = minimalDistance, 
    fun_R = meanFun),
  type = "costFun.optimiseSD"
  )[[1]]

# big example: slow
optSD_manual1 = optimiseSD_manual(simulations = radioactivePlumes_local, 
                             costFun = spatialSpreadMinDist, 
                             costMap = spatialSpreadMinDist,
                             locationsAll = locAll1, 
                             locationsFix = locFix1, 
                             locationsInitial = locInit1)

spatialSpreadMinDist1 = spatialSpreadMinDist(simulations = radioactivePlumes_loc, locations = 1)
optSD_manual2 = optimiseSD_manual(simulations = radioactivePlumes_loc, 
                                  costFun = spatialSpreadMinDist, 
                                  costMap = spatialSpreadMinDist,
                                  locationsAll = 1991:2001,#locAll1, 
                                  locationsFix = 1:1990,#locFix1, 
                                  locationsInitial = locInit1)
optSD_manual2 = optimiseSD_manual(simulations = radioactivePlumes_loc, 
                                  costFun = spatialSpreadMinDist, 
                                  costMap = spatialSpreadMinDist, 
                                  locationsFix = locFix1, 
                                  locationsInitial = locInit1)
# - - - - - - - - - SIndexDF - - - - - - -
data(SIndexDF)
simulationsIndex = subset(radioactivePlumes_local, locations = c(220:222), plumes = 11:15)
simulationsIndex@locations = SIndexDF
plot(simulationsIndex)
spplot(simulationsIndex@locations)
## layer: detection 
simulationsIndex@values$detectable = calc(
  simulationsIndex@values$maxdose,
  fun = function(x){x >= 1e-2}) # threshold chosen to get interesting result
## weight
simulationsIndex@plumes$totalDose = 
  summaryPlumes(simulationsIndex, fun = sum, values = "finaldose")[[2]]

singleDet1 = singleDetection(simulations = simulationsIndex, locations = 1)
#setBreakpoint("optimiseSD_manual.R#146")
optSD_manual3 = optimiseSD_manual(simulations = simulationsIndex, 
                                  costFun = singleDetection, 
                                  costMap = singleDetection,
                                  valuesPlot = "a",
                                  col = terrain.colors)


# - - - - - - - SPixelsDF - - - - - -
data(SPixelsDF)
simulationsPixels = subset(radioactivePlumes_local, 
                           locations = c(199:201, 220:222, 241:243), plumes = 21:25)
simulationsPixels@locations = SPixelsDF
plot(simulationsPixels)
spplot(simulationsPixels@locations)

## layer: detection 
simulationsPixels@values$detectable = calc(
  simulationsPixels@values$maxdose,
  fun = function(x){x >= 1e-2}) # threshold chosen to get interesting result
plot(simulationsPixels, zcol = "detectable")
## weight
simulationsPixels@plumes$totalDose = 
  summaryPlumes(simulationsPixels, fun = sum, values = "finaldose")[[2]]

singleDet_Pixels89 = singleDetection(simulations = simulationsPixels,locations = c(8,9))
optSD_manual4 = optimiseSD_manual(simulations = simulationsPixels, 
                                  locationsFix = 9,
                                  costFun = singleDetection, 
                                  costMap = singleDetection,
                                  col = heat.colors)

# - - - - - - - SPointsDF - - - - - - - - - - -
data(SPointsDF)
simulationsPoints = subset(radioactivePlumes_local, 
                           locations = c(199:201, 220:222, 241:243), plumes = 21:25)
simulationsPoints@locations = SPointsDF
plot(simulationsPoints)
spplot(simulationsPoints@locations)
simulationsPoints@values$detectable = calc(
  simulationsPoints@values$maxdose,
  fun = function(x){x >= 1e-2}) # threshold chosen to get interesting result
plot(simulationsPoints, zcol = "detectable")
simulationsPoints@plumes$totalDose = 
  summaryPlumes(simulationsPoints, fun = sum, values = "finaldose")[[2]]

singleDet_Points89 = singleDetection(simulations = simulationsPoints, locations = c(8,9))
optSD_manual5 = optimiseSD_manual(simulations = simulationsPoints, 
                                  locationsFix = 3,
                                  costFun = singleDetection, 
                                  costMap = singleDetection,
                                  col = topo.colors
                                  )



# - - - - - - - - - SPolygonsDF - - - - - - - - - - 
data(SPolygonsDF)
simulationsPolygons = subset(radioactivePlumes_local, 
                           locations = c(200, 221, 242), plumes = 31:35)
simulationsPolygons@locations = SPolygonsDF
plot(simulationsPolygons)
spplot(simulationsPolygons@locations)
simulationsPolygons@values$detectable = calc(
  simulationsPolygons@values$maxdose,
  fun = function(x){x >= 1e-3}) # threshold chosen to get interesting result
plot(simulationsPolygons, zcol = "detectable")
simulationsPolygons@plumes$totalDose = 
  summaryPlumes(simulationsPolygons, fun = sum, values = "finaldose")[[2]]

singleDet_Polygons2 = singleDetection(simulations = simulationsPolygons, locations = 2)
optSD_manual7 = optimiseSD_manual(simulations = simulationsPolygons, 
                                  costFun = singleDetection, 
                                  costMap = singleDetection,
                                  col = heat.colors)


# - - - - - - - SLinesDF - - - - - - - - - -
data(SLinesDF)
simulationsLines = subset(radioactivePlumes_local, 
                             locations = c(200, 221), plumes = 41:45)
simulationsLines@locations = SLinesDF
plot(simulationsLines)
spplot(simulationsLines@locations)
simulationsLines@values$detectable = calc(
  simulationsLines@values$maxdose,
  fun = function(x){x >= 1e-3}) # threshold chosen to get interesting result
plot(simulationsLines, zcol = "detectable")
simulationsLines@plumes$totalDose = 
  summaryPlumes(simulationsLines, fun = sum, values = "finaldose")[[2]]

singleDet_Lines2 = singleDetection(simulations = simulationsLines, locations = 2)
#setBreakpoint("optimiseSD_manual.R#146")
optSD_manual8 = optimiseSD_manual(simulations = simulationsLines, 
                                  costFun = singleDetection, 
                                  costMap = singleDetection,
                                  col = topo.colors)

# - - - - - - - - via optimiseSD - - - -- - - - - -  
optSDman_singleDet = replaceDefault(
  optimiseSD_manual,
  newDefaults = list(
    costMap = singleDetection,
    col = topo.colors
    ),
  type = "optimisationFun.optimiseSD")[[1]]

optSD_manual9 = optimiseSD(
  simulations = simulationsLines, 
  costFun = singleDetection,
  optimisationFun = optSDman_singleDet)
