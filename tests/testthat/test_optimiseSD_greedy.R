################################################################
# test optimiseSD_greedy                                       #
################################################################
# deleteSensor
# addSensor
# determineNextStep
data(SimulationsSmall)
locDel1 = sample.int(nLocations(SimulationsSmall), 2)
locKeep1 = sample(setdiff(1:nLocations(SimulationsSmall), locDel1), 2)
locAll1 = c(sample(setdiff(1:nLocations(SimulationsSmall), c(locDel1, locKeep1)), 4), locDel1)
data(radioactivePlumes_local)
locDel2 = sample.int(nLocations(radioactivePlumes_local), 5)
locKeep2 = sample(setdiff(1:nLocations(radioactivePlumes_local), locDel2), 100)
locAll2 = c(sample(setdiff(1:nLocations(radioactivePlumes_local), c(locDel2, locKeep2)), 10), locDel2)
data(radioactivePlumes_area)
locDel3 = sample.int(nLocations(radioactivePlumes_area), 5)
locKeep3 = sample(setdiff(1:nLocations(radioactivePlumes_area), locDel3), 100)
locAll3 = c(sample(setdiff(1:nLocations(radioactivePlumes_area), c(locDel3, locKeep3)), 10), locDel3)

meanFun = function(x){mean(x, na.rm = TRUE)}
minDist = replaceDefault(
  spatialSpread, newDefaults = list(
    fun = minimalDistance,
    fun_R = meanFun
  ), type = "costFun.optimiseSD"
)[["fun"]] 
# ---------------------------------- deleteSensor ------------------------ #
del1 = deleteSensor(simulations = SimulationsSmall, 
             costFun = minDist, 
             locationsDeletable = locDel1, 
             locationsKeep = locKeep1)

costWithout1 = numeric(length(locDel1))
for (i in seq(along = locDel1)){
  costWithout1[i] = minDist(SimulationsSmall, c(locDel1[-i], locKeep1))[[1]]
}
expect_equal(
  del1[[3]],
  costWithout1)
expect_equal(
  del1[[2]],
  which(costWithout1 == min(costWithout1))
  )
expect_true(
  is.element(setdiff(locDel1, del1[[1]]), locDel1[which(costWithout1 == min(costWithout1))])
  )

#
del2 = deleteSensor(simulations = radioactivePlumes_local, 
                    costFun = minDist, 
                    locationsDeletable = locDel2, 
                    locationsKeep = locKeep2)

costWithout2 = numeric(length(locDel2))
for (i in seq(along = locDel2)){
  costWithout2[i] = minDist(radioactivePlumes_local, c(locDel2[-i], locKeep2))[[1]]
}
expect_equal(
  del2[[3]],
  costWithout2)
expect_equal(
  del2[[2]],
  which(costWithout2 == min(costWithout2))
)
expect_true(
  is.element(setdiff(locDel2, del2[[1]]), locDel2[which(costWithout2 == min(costWithout2))])
)

# ---------------- addSensor ----------------------------- #
add_delinE_kl1 = addSensor(
  costFun = delineationError_kl1,                                                            
  simulations = radioactivePlumes_local,
  locationsAddable = locDel2,        # currently empty potential sensor locations; one of them shall be added
  locationsCurrent = locKeep2         # current sampling design
)

costWith1 = numeric(length(locDel2))
for (i in seq(along = locDel2)){
  costWith1[i] = delineationError_kl1(radioactivePlumes_local, c(locDel2[i], locKeep2))[[1]]
}
expect_equal(
  add_delinE_kl1[[3]],
  costWith1)
expect_equal(
  add_delinE_kl1[[2]],
  which(costWith1 == min(costWith1))
)
expect_true(
  is.element(add_delinE_kl1[[1]], locDel2[which(costWith1 == min(costWith1))])
)

# ------------------------ determineNextStep ------------------------ #
# - - - - -  kindAim = "number" - - - - - 
# stop because aim reached (aimNumber = length(locCurrent))
next_stopIt = determineNextStep(
  locationsCurrent = locDel1,                      # indices of current sensors
  locationsFix = locKeep1,                      # indices of current sensors that are fix
  locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
  kindAim = "number",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
  aim = 4,                                        # desired number of sensors or desired cost
  maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria                            
  l = 10,                                            # iteration step
  costFun = minDist,                                # to be forwarded to computeCost
  simulations = SimulationsSmall
)
expect_equal(next_stopIt[[1]], "finish")

# del, because aim not reached (aimNumber < length(locCurrent))
next_stopIt = determineNextStep(
  locationsCurrent = locDel1,                      # indices of current sensors
  locationsFix = locKeep1,                      # indices of current sensors that are fix
  locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
  kindAim = "number",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
  aim = 3,                                        # desired number of sensors or desired cost
  maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria                            
  l = 10,                                            # iteration step
  costFun = minDist,                                # to be forwarded to computeCost
  simulations = SimulationsSmall
)
expect_equal(next_stopIt[[1]], "del")

# add, because aim not reached (aimNumber > length(locCurrent))
next_stopIt = determineNextStep(
  locationsCurrent = locDel1,                      # indices of current sensors
  locationsFix = locKeep1,                      # indices of current sensors that are fix
  locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
  kindAim = "number",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
  aim = 5,                                        # desired number of sensors or desired cost
  maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria                            
  l = 10,                                            # iteration step
  costFun = minDist,                                # to be forwarded to computeCost
  simulations = SimulationsSmall
)
expect_equal(next_stopIt[[1]], "add")

# - - - - -  kindAim = "cost" - - - - - 
# stop because aim reached (aimCost = cost(locCurrent))
next_stopIt = determineNextStep(
  locationsCurrent = locDel1,                      # indices of current sensors
  locationsFix = locKeep1,                      # indices of current sensors that are fix
  locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
  kindAim = "cost",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
  aim = minDist(simulations = SimulationsSmall,
                locations = locDel1)[[1]],                                        # desired number of sensors or desired cost
  maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria                            
  l = 10,                                            # iteration step
  costFun = minDist,                                # to be forwarded to computeCost
  simulations = SimulationsSmall
)
expect_equal(next_stopIt[[1]], "finish")

# del because aim not reached (aimCost > cost(locCurrent))
next_stopIt = determineNextStep(
  locationsCurrent = locDel1,                      
  locationsFix = locKeep1,                      
  locationsAll = locAll1,            
  kindAim = "cost",                               
  aim = minDist(simulations = SimulationsSmall,
                locations = locDel1)[[1]] * 2,                                        
  maxIterations = 10,                       
  l = 10,                                            # iteration step
  costFun = minDist,                                # to be forwarded to computeCost
  simulations = SimulationsSmall
)
expect_equal(next_stopIt[[1]], "del")

# add, because aim not reached (aimCost < cost(locCurrent))
next_stopIt = determineNextStep(
  locationsCurrent = locDel1,                      # indices of current sensors
  locationsFix = locKeep1,                      # indices of current sensors that are fix
  locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
  kindAim = "cost",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
  aim = minDist(simulations = SimulationsSmall,
                locations = locDel1)[[1]] * 1/2,                                        # desired number of sensors or desired cost
  maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria                            
  l = 10,                                            # iteration step
  costFun = minDist,                                # to be forwarded to computeCost
  simulations = SimulationsSmall
)
expect_equal(next_stopIt[[1]], "add")

# stop
## stop because iterations full (l > maxIterations)
expect_warning(
  next_stopIt <- determineNextStep(
    locationsCurrent = locDel1,                      # indices of current sensors
    locationsFix = locKeep1,                      # indices of current sensors that are fix
    locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
    kindAim = "number",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
    aim = 4,                                        # desired number of sensors or desired cost
    maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria                            
    l = 12,                                            # iteration step
    costFun = minDist,                                # to be forwarded to computeCost
    simulations = SimulationsSmall
  )  
)
expect_equal(next_stopIt[[1]], "stop")

## stop because all sensors used and cost not reached
expect_warning(
  next_stopIt <- determineNextStep(
    locationsCurrent = locAll1,                      # indices of current sensors
    locationsFix = locKeep1,                      # indices of current sensors that are fix
    locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
    kindAim = "cost",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
    aim = 0,                                        # desired number of sensors or desired cost
    maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria                            
    l = 10,                                            # iteration step
    costFun = minDist,                                # to be forwarded to computeCost
    simulations = SimulationsSmall
  )  
)
expect_equal(next_stopIt[[1]], "stop")

## stop because all sensors included and number not reached
expect_warning(
  next_stopIt <- determineNextStep(
    locationsCurrent = 1:9,                      
    locationsFix = locKeep1,                      
    locationsAll = locAll1,           
    kindAim = "number",                               
    aim = 10,                                        
    maxIterations = 10,                      
    l = 10,                                            
    costFun = minDist,                                
    simulations = SimulationsSmall
  )  
)
expect_equal(next_stopIt[[1]], "stop")

## stop because only fix sensors used and cost still below aimCost
expect_warning(
  next_stopIt <- determineNextStep(
    locationsCurrent = locDel1,                      # indices of current sensors
    locationsFix = locDel1,                      # indices of current sensors that are fix
    locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
    kindAim = "cost",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
    aim = minDist(
      simulations = SimulationsSmall,
      locations = locDel1[-1])[[1]],                                        # desired number of sensors or desired cost
    maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria                            
    l = 10,                                            # iteration step
    costFun = minDist,                                # to be forwarded to computeCost
    simulations = SimulationsSmall
  )
)
expect_equal(next_stopIt[[1]], "stop")


## stop because only fix sensors used and number still above aimNumver
expect_warning(
  next_stopIt <- determineNextStep(
    locationsCurrent = locDel1,                      # indices of current sensors
    locationsFix = locDel1,                      # indices of current sensors that are fix
    locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
    kindAim = "number",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
    aim = 3,                                        # desired number of sensors or desired cost
    maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria                            
    l = 10,                                            # iteration step
    costFun = minDist,                                # to be forwarded to computeCost
    simulations = SimulationsSmall
  )
)
expect_equal(next_stopIt[[1]], "stop")


# ---------------------- optimiseSD_greedy ---------------------------- #
# aimNumber (exact), add, no swap
optGreedy1 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimNumber = 8,
  maxIterations = 20,             
  swap = FALSE                      
)
expect_equal(
  optGreedy1$evalSDs$number,
  6:8
)
expect_equal(
  optGreedy1$finalSDwhich,
  3
)
# aimNumber (exact), delete, no swap
optGreedy2 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimNumber = 4,
  maxIterations = 20,             
  swap = FALSE                      
)
expect_equal(
  optGreedy2$evalSDs$number,
  6:4
)
expect_equal(
  optGreedy2$finalSDwhich,
  3
)
# aimCost (to undershoot), add, no swap
optGreedy3 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimCost = 0.17,
  maxIterations = 20,             
  swap = FALSE                      
)
expect_equal(optGreedy1, optGreedy3)
# aimCost (to undershoot), del, no swap
optGreedy4 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimCost = 1.3,
  maxIterations = 20,             
  swap = FALSE                      
)
expect_equal(
  optGreedy4$evalSDs$number,
  6:3
)
expect_equal(
  optGreedy4$finalSDwhich,
  3
)
# aimNumber(exact), add, swap
optGreedy5 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimNumber = 8,
  maxIterations = 20,             
  swap = TRUE                      
)
expect_equal(
  optGreedy5$evalSDs$number,
  c(6:8, 7)
)
expect_equal(
  optGreedy5$finalSDwhich,
  3
)
# aimNumber(exact), delete, swap
optGreedy6 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimNumber = 4,
  maxIterations = 20,             
  swap = TRUE                      
)
expect_equal(
  optGreedy6$evalSDs$number,
  c(6:4, 5, 4, 5) # why two tries?
)
expect_equal(
  optGreedy6$finalSDwhich,
  c(3, 5)
)
# aimCost (to undershoot), add, swap
optGreedy7 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimCost = 0.17,
  maxIterations = 20,             
  swap = TRUE                      
)
expect_equal(
  optGreedy7$evalSDs$number,
  c(6:8,7) 
)
expect_equal(
  optGreedy7$finalSDwhich,
  3
)
# aimCost (to undershoot), del, swap
optGreedy8 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimCost = 1.3,
  maxIterations = 20,             
  swap = TRUE                      
)
expect_equal(
  optGreedy8$evalSDs$number,
  c(6:3,4,3,4) # why two tries?
)
expect_equal(
  optGreedy8$finalSDwhich,
  c(3,5,7)
)
# various reasons to stop (with/without swap)
# cost, add, maxIterations, no swap
optGreedy9 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimCost = 0.17,
  maxIterations = 2,             
  swap = FALSE                      
)
expect_equal(
  optGreedy9$evalSDs$number,
  c(6:8) 
)
expect_equal(
  optGreedy9$finalSDwhich,
  3
)
# cost, del, not reachable, no swap
optGreedy10 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimCost = 0,
  maxIterations = 20,             
  swap = FALSE                      
)
expect_equal(
  optGreedy10$evalSDs$number,
  c(6:8) 
)
expect_equal(
  optGreedy10$finalSDwhich,
  integer(0)
)
# number, add, not reachable, no swap
optGreedy11 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimNumber = 9,
  maxIterations = 20,             
  swap = FALSE                      
)
expect_equal(
  optGreedy11$evalSDs$number,
  c(6:8) 
)
expect_equal(
  optGreedy11$finalSDwhich,
  3
)
# number, del, maxIterations, swap
optGreedy12 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimNumber = 4,
  maxIterations = 1,             
  swap = TRUE                      
)
expect_equal(
  optGreedy12$evalSDs$number,
  c(6:5) 
)
expect_equal(
  optGreedy12$finalSDwhich,
  integer(0)
)
# cost, del, not reachable, no swap
optGreedy13 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimCost = 4,
  maxIterations = 20,             
  swap = FALSE                      
)
expect_equal(
  optGreedy13$evalSDs$number,
  c(6:2) 
)
expect_equal(
  optGreedy13$finalSDwhich,
  5
)
# number, add, not reachable, swap
optGreedy14 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimNumber = 9,
  maxIterations = 20,             
  swap = TRUE                      
)
expect_equal(
  optGreedy14$evalSDs$number,
  c(6:8) 
)
expect_equal(
  optGreedy14$finalSDwhich,
  3
)
# number (ignored) and cost, swap
optGreedy15 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimCost = 1.3,
  aimNumber = 9,
  maxIterations = 20,             
  swap = TRUE                      
)
expect_equal(
  optGreedy15,
  optGreedy8
)


# big dataset
optGreedy_l1 = optimiseSD_greedy(
  simulations = radioactivePlumes_local,
  costFun = minDist,                        
  locationsAll = locAll2,                    
  locationsFix = locKeep2,                   
  locationsInitial = locDel2,                
  aimNumber = 101,
  maxIterations = 20,             
  swap = TRUE                      
)
optGreedy_l2 = optimiseSD_greedy(
  simulations = radioactivePlumes_local,
  costFun = minDist,                        
  locationsAll = locAll2,                    
  locationsFix = locKeep2,                   
  locationsInitial = locDel2,                
  aimNumber = 109,
  maxIterations = 20,             
  swap = TRUE                      
)
# saving
# file.remove("data/optGr8.Rdata")
# filesBefore = list.files("data")
# optGreedy8_ = optimiseSD_greedy(
#   simulations = SimulationsSmall,
#   costFun = minDist,                        
#   locationsAll = locAll1,                    
#   locationsFix = locKeep1,                   
#   locationsInitial = locDel1,                
#   aimCost = 1.3,
#   maxIterations = 20,             
#   swap = TRUE,
#   nameSave = "data/optGr8"
# )
# filesAfter = list.files("data")
# setdiff(filesAfter, filesBefore) == "optGr8.Rdata"
# load("data/optGr8.Rdata")
# expect_equal(
#   SDs,
#   optGreedy8_$SDs
# )
# file.remove("data/optGr8.Rdata")

# plotting: was parameter to be forwarded to cost function each time it is called; else no functionality

# optimiseSD_greedy called by optimiseSD

optGreedy8 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,                        
  locationsAll = locAll1,                    
  locationsFix = locKeep1,                   
  locationsInitial = locDel1,                
  aimCost = 1.3,
  maxIterations = 20,             
  swap = TRUE                      
)

optimSD_greedy = replaceDefault(
  optimiseSD_greedy,  newDefaults = list(
    maxIterations = 20,             
    swap = TRUE    
    ),
  type = "optimisationFun.optimiseSD")[[1]] 

opt8 = optimiseSD(
  simulations = SimulationsSmall,                                      
  costFun = minDist,
  locationsAll = locAll1,        
  locationsFix = locKeep1,
  locationsInitial = locDel1,                   
  aimCost = 1.3,
  optimisationFun = optimSD_greedy,
  nameSave = "data/optimiseSD8"
)
expect_equal(opt8, optGreedy8)

optimSD_greedy2 = replaceDefault(
  optimiseSD_greedy,  newDefaults = list(
    maxIterations = 20,             
    swap = FALSE   
  ),
  type = "optimisationFun.optimiseSD")[[1]] 

opt4 = optimiseSD(
  simulations = SimulationsSmall,                                      
  costFun = minDist,
  locationsAll = locAll1,        
  locationsFix = locKeep1,
  locationsInitial = locDel1,                   
  aimCost = 1.3,
  optimisationFun = optimSD_greedy2,
  nameSave = "data/optimiseSD8"
)
expect_equal(opt4, optGreedy4)
