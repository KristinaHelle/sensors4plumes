######################################################################
# test plotSD                                                        #
######################################################################

# ssa
## based on example from optimiseSD_ssa manual
if (prepare = TRUE){
  # the function is to be used inside of optimiseSD
  # change parameters 
  optimSD_ssa1 = replaceDefault(
    optimiseSD_ssa,  newDefaults = list(
      start_acc_vG = 0.1,
      aimCost = 0, 
      verbatim = TRUE,
      maxIterations = 3000,
      maxStableIterations = 500,
      maxIterationsJumpBack = 200   
    ),
    type = "optimisationFun.optimiseSD")[[1]] 
  
  # define possible, fix, and initial sensors
  data(radioactivePlumes_area)
  locDel3 = sample.int(nLocations(radioactivePlumes_area), 5)
  locKeep3 = sample(setdiff(1:nLocations(radioactivePlumes_area), locDel3), 100)
  locAll3 = c(sample(setdiff(1:nLocations(radioactivePlumes_area), 
                             c(locDel3, locKeep3)), 10), locDel3)
  
  # prepare data
  threshold = 1e-7
  radioactivePlumes_area@values$detectable = calc(
    radioactivePlumes_area@values$maxdose,
    fun = function(x){x >= threshold})
  radioactivePlumes_area@plumes$nDetectable = 
    summaryPlumes(radioactivePlumes_area, fun = sum, values = "detectable")[[2]]  
  radioactivePlumes_area@plumes$totalDose = 
    summaryPlumes(radioactivePlumes_area, fun = sum, values = "finaldose")[[2]]
  radioactivePlumes_area@plumes$earliestDetection = 
    summaryPlumes(radioactivePlumes_area, fun = min, values = "time", na.rm = TRUE)[[2]]
  
  costInitial1 = multipleDetection(simulations = radioactivePlumes_area, 
                                   locations = c(locKeep3, locDel3))
  
  # run optimisation
  optSSA1 = optimiseSD(
    simulations = radioactivePlumes_area,
    costFun = singleDetection,
    locationsAll = setdiff(1:nLocations(radioactivePlumes_area), c(locKeep3, locAll3)),                       
    locationsFix = locKeep3,
    locationsInitial = locDel3,
    aimCost = 0.05 * costInitial1[[1]], 
    aimNumber = length(locDel3),
    optimisationFun = optimSD_ssa1
  )   
}
# one SD, no costMap
plotSSA1_0 = plotSD(
  simulations = radioactivePlumes_area,
  SD = optSSA1$SD,
  locationsAll = setdiff(1:nLocations(radioactivePlumes_area), c(locKeep3, locAll3)),                       
  locationsFix = locKeep3,
  locationsInitial = locDel3
)
plotSSA1_0
# one SD, no costMap, customised plotting
plotSSA1_1 = plotSD(
    simulations = radioactivePlumes_area,
    SD = optSSA1$SD,
    locationsAll = setdiff(1:nLocations(radioactivePlumes_area), c(locKeep3, locAll3)),                       
    locationsFix = locKeep3,
    locationsInitial = locDel3,
    pch = c(4, 10, 12),
    col = c("green", "white", "lightgreen", "cyan"),
    col.regions = heat.colors,
    colorkey = FALSE
  )
plotSSA1_1

# one SD, costMap
plotSSA1_2 = plotSD(
  simulations = radioactivePlumes_area,
  SD = optSSA1$SD,
  locationsAll = setdiff(1:nLocations(radioactivePlumes_area), c(locKeep3, locAll3)),                       
  locationsFix = locKeep3,
  locationsInitial = locDel3,
  costMap = singleDetection,
  cuts = 3
)
plotSSA1_2


## prepare costMap
meanFun = function(x){mean(x, na.rm = TRUE)}
spatialSpread_minDist = replaceDefault(
  spatialSpread,
  newDefaults = list(
    weightByArea = TRUE,
    fun = minimalDistance,
    fun_R = meanFun
    )
)[[1]]
# SDs matrix
SDs = matrix(sample.int(2500, 40), nrow = 8)
# separate maps
plotSSA1_2 = plotSD(
  simulations = radioactivePlumes_area,
  SD = SDs[1:6,],                      
  locationsFix = SDs[7,],
  locationsInitial = SDs[8,],
  costMap = spatialSpread_minDist,
  at = c(0,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000),
  mainCost = FALSE
)
c(plotSSA1_2[[1]], plotSSA1_2[[2]], plotSSA1_2[[3]], plotSSA1_2[[4]], plotSSA1_2[[5]], plotSSA1_2[[6]],
  layout = c(2,3)) # one key, no main

plotSSA1_2a = plotSD(
  simulations = radioactivePlumes_area,
  SD = SDs[1:6,],                      
  locationsFix = SDs[7,],
  locationsInitial = SDs[8,],
  costMap = spatialSpread_minDist
)
# par(mfrow = c(2,3)) # this line would have no effect, generates single maps
 for (i in 1:6){
   plot(plotSSA1_2a[[i]])
 }

# single map
plotSSA1_3 = plotSD(
  simulations = radioactivePlumes_area,
  SD = SDs[1:6,],                      
  locationsFix = SDs[7,],
  locationsInitial = SDs[8,],
  costMap = spatialSpread_minDist,
  allIn1Plot = 4,
  pch = c(1,20,4)
)
plotSSA1_3 

# single map, customise style for SDs
plotSSA1_3 = plotSD(
  simulations = radioactivePlumes_area,
  SD = SDs[1:6,],                      
  locationsFix = SDs[7,],
  costMap = spatialSpread_minDist,
  allIn1Plot = 4,
  pch = c(1,20,19),
  pch.SDs = 1:6,
  col.SDs = 1:6,
  cex.SDs = rep(1, times = 6),
  col.regions = grey.colors,
  pointsKey = FALSE
)
plotSSA1_3 

# SDs list
SDs2 = list(
  sample.int(2500, 5),
  sample.int(2500, 10),
  sample.int(2500, 15),
  sample.int(2500, 20)
  )
plotSSA1_4 = plotSD(
  simulations = radioactivePlumes_area,
  SD = SDs2,                      
  costMap = spatialSpread_minDist,
  allIn1Plot = 2
)
plotSSA1_4

