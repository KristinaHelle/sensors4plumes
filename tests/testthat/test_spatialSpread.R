################################################
# test spread cost functions                   #
################################################


# test with big number of locations and sensors
# test without locations, wrong locations, locations of length 0, 1, many

testthat("spatialSpread", {
  data(radioactivePlumes_local)
   
  meanFun = function(x){mean(x, na.rm = TRUE)}
  

  locations1 = sample.int(nLocations(radioactivePlumes_local), 7000)
  locations2 = sample.int(nLocations(radioactivePlumes_local), 1)
  locations3 = sample.int(nLocations(radioactivePlumes_local), 200)
  # minimalDistance
  ## many locations  
  spatialSpread1 = spatialSpread(
    simulations = radioactivePlumes_local,
    locations = c(NA, 0, 8000, locations1, rev(locations1)), 
    weightByArea = TRUE,
    fun = minimalDistance,
    fun_R = meanFun
  )
  expect_true(
    all(spatialSpread1[["costLocations"]][locations1] == 0)
  )
  expect_true(
    all(spatialSpread1[["costLocations"]][-locations1] > 0)
  )
  expect_equivalent(
    spatialSpread1[["costLocations"]],
    get.knnx(coordinates(radioactivePlumes_local@locations)[locations1,],
             coordinates(radioactivePlumes_local@locations),
             algorithm = "kd_tree")$nn.dist[,1]  
    )
  
  expect_equivalent(
    spatialSpread1[["cost"]],
    mean(spatialSpread1[["costLocations"]] * 
           areaSDF(radioactivePlumes_local@locations))  
  )
  
  ## 1 location
  spatialSpread2 = spatialSpread(
    simulations = radioactivePlumes_local,
    locations = locations2, 
    weightByArea = TRUE,
    fun = minimalDistance,
    fun_R = meanFun
  )
  expect_equivalent(
    spatialSpread2[["costLocations"]],
    get.knnx(coordinates(radioactivePlumes_local@locations)[locations2,,drop = FALSE],
             coordinates(radioactivePlumes_local@locations), k = 1,
             algorithm = "kd_tree")$nn.dist[,1]  
  )
  ## 0 location
  spatialSpread3 = spatialSpread(
    simulations = radioactivePlumes_local,
    locations = integer(0), 
    weightByArea = TRUE,
    fun = minimalDistance,
    fun_R = meanFun
  )
  expect_true(
    all(spatialSpread3[["costLocations"]] == Inf)
  )

  # krigingVariance
  samplePlumes = sort(sample.int(nPlumes(radioactivePlumes_local), 100))
  medianVariogram = fitMedianVariogram(simulations = radioactivePlumes_local, 
                                       plumes = samplePlumes, 
                                       locations = locations3, # using all data takes too long
                                       values = 1)
  krigingVarianceMedian = replaceDefault(krigingVariance, 
                                         newDefaults = list(model = medianVariogram))[["fun"]]
  # 1 locations
  spatialSpread4 = spatialSpread(
    simulations = radioactivePlumes_local,
    locations = c(NA, 0, 8000, locations2, rev(locations2)), 
    weightByArea = TRUE,
    fun = krigingVarianceMedian,
    fun_R = meanFun
  )

  # 0 locations
  spatialSpread5 = spatialSpread(
    simulations = radioactivePlumes_local,
    locations = integer(0), 
    weightByArea = TRUE,
    fun = krigingVarianceMedian,
    fun_R = meanFun
  )
  expect_true(
    all(spatialSpread5[["costLocations"]] == Inf)
  )
  spatialSpread6 = spatialSpread(
    simulations = radioactivePlumes_local,
    locations = locations2, 
    weightByArea = TRUE,
    fun = krigingVarianceMedian,
    fun_R = meanFun
  )
  
  spatialSpread7 = spatialSpread(
    simulations = radioactivePlumes_local,
    locations = locations3, 
    weightByArea = TRUE,
    fun = krigingVarianceMedian,
    fun_R = meanFun
  )
  expect_true( #numeric error, should be 0
    all(abs(spatialSpread7[["costLocations"]][locations3]) < 1e-17)
  )
  expect_true(
    all(spatialSpread7[["costLocations"]][-locations3] > 0)
  )
  
  asPoints = as(radioactivePlumes_local@locations, "SpatialPointsDataFrame")
  kV7 = krige(index ~ 1, asPoints[locations3,], asPoints, model =  medianVariogram)
  
  expect_equal(
    kV7$var1.var,
    spatialSpread7[["costLocations"]]
  )
})
