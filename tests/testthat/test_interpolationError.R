###########################################################################
# test errorInterpolation                                                 #
###########################################################################

# test the single error functions
test_that("errorInterpolationFunctions", {
  X = c(0.9 ,3.5)
  expect_equal(absError(x = X), 
               abs(diff(X)))
  expect_equal(delineationError(x = X, threshold = 1)[3], 1)
  expect_equal(delineationError(x = rev(X), threshold = 1)[3], 5)
  expect_equal(delineationError(x = X, threshold = 0.5)[3], 0)
})


# test errorInterpolation altogether
test_that("interpolationError", {
  data(radioactivePlumes_area)
  ## preparation
  samplePlumes = sort(sample.int(nPlumes(radioactivePlumes_area), 100))
  medianVariogram = fitMedianVariogram(simulations = radioactivePlumes_area, 
                                       plumes = samplePlumes, 
                                       values = 1)
  krige0var = replaceDefault(krige0, newDefaults = list(
    formula = z ~ 1, model = medianVariogram, beta = NA, ... = NA))[[1]]
  idw0z = replaceDefault(idw0, newDefaults = list(
    formula = z ~ 1))[[1]]
  sampleLocations250 = sample.int(nLocations(radioactivePlumes_area), 250)
  fun_Rpl_mean3 = function(x, nout = 3){ 
    c(mean(x[,1], na.rm = TRUE),
      mean(x[,2], na.rm = TRUE),
      mean(x[,3], na.rm = TRUE))
  }
  fun_Rpl_mean1 = function(x, nout = 1){ 
    c(mean(x[,1], na.rm = TRUE))
  }
  ## interpolated / error in memory
  ## fun_Rpl given
  expect_warning(
    interpolationError1 <- interpolationError(
      simulations = radioactivePlumes_area,
      locations = c(NA, NA, 0, 0, 0, sampleLocations250, rev(sampleLocations250)),
      values = 2,
      fun_interpolation = krige0var,
      fun_error = delineationError,
      fun_l = delineationErrorMap,
      fun_Rpl = fun_Rpl_mean3,
      fun_Rpl_cellStats = "mean",
      tmpfile = "interpolationError1",
      overwrite = TRUE,
      chunksize = 1e+7
    )    
  )

  ## fun_Rpl invalid
  expect_warning(
    interpolationError2 <- interpolationError(
      simulations = radioactivePlumes_area,
      locations = c(NA, NA, 0, 0, 0, sampleLocations250, rev(sampleLocations250)),
      values = 2,
      fun_interpolation = krige0var,
      fun_error = delineationError,
      fun_Rpl = "fun_Rpl_mean3",
      fun_Rpl_cellStats = "mean",
      tmpfile = "interpolationError2",
      overwrite = TRUE,
      chunksize = 1e+7
    )    
  )
  
  expect_warning(# other interpolation and error function
    interpolationError3 <- interpolationError(
      simulations = radioactivePlumes_area,
      locations = sampleLocations250,
      values = 2,
      fun_interpolation = idw0z,
      fun_error = absError,
      fun_l = absErrorMap,
      fun_Rpl = fun_Rpl_mean1,
      fun_Rpl_cellStats = "mean",
      tmpfile = "interpolationError3",
      overwrite = TRUE,
      chunksize = 1e+7
    )    
  )
  ## for comparison
  interpolation1 = interpolate(
    simulations = radioactivePlumes_area,
    locations = c(NA, NA, 0, 0, 0, sampleLocations250, rev(sampleLocations250)),
    values = 2,
    fun_interpolation = krige0var
    )
  expect_equal(# "cost"
    interpolationError1[["cost"]],
    fun_Rpl_mean3(getValues(interpolationError1[["error_locationsplumes"]]))
  )
#  expect_equal(#"cost_cellStats" (fails because needs numeric tolerance ~ 1e-16)
#    interpolationError1[["cost_cellStats"]],
#    interpolationError1[["cost"]]
#  )
  
  i = sample.int(nLocations(radioactivePlumes_area), 1)
  j = sample.int(nPlumes(radioactivePlumes_area), 1)  
  expect_equivalent(# "error_locationsplumes"
    delineationError(
      x = c(interpolationError1[["interpolated"]][i,j],
            radioactivePlumes_area@values[[2]][i,j])),
    as.numeric(interpolationError1[["error_locationsplumes"]][i,j,])
  )
  expect_equal(# "interpolated"
    interpolationError1[["interpolated"]],
    interpolation1)

  expect_equal(# 'costMap'
    interpolationError1[["costLocations"]][i],
    mean(interpolationError1[["error_locationsplumes"]][i,][,1] == 0 &
          interpolationError1[["error_locationsplumes"]][i,][,2] == 1) + 
      5 * 
      mean(interpolationError1[["error_locationsplumes"]][i,][,1] == 1 &
            interpolationError1[["error_locationsplumes"]][i,][,2] == 0)
  )
  expect_equal(# 'costMap'
    interpolationError1[["costLocations"]][i],
    mean(interpolationError1[["error_locationsplumes"]][i,][,3])
  )

  # if no valid fun_Rpl, 'cost' is 'cost_cellStats'
  expect_equivalent(
    interpolationError1[2:4],
    interpolationError2
  )
  expect_equal(
    names(interpolationError1)[c(-2, -5)],
    names(interpolationError2)
  )

  ## interpolated / error not in memory
  sampleLocations2000 = sample.int(nLocations(radioactivePlumes_area), 2000)

  expect_warning(
    interpolationError4 <- interpolationError(
      simulations = radioactivePlumes_area,
      locations = c(NA, NA, 0, 0, 0, sampleLocations2000, rev(sampleLocations2000)),
      values = 2,
      fun_interpolation = krige0var,
      fun_error = delineationError,
      fun_Rpl = fun_Rpl_mean3,
      fun_Rpl_cellStats = "mean",
      tmpfile = "interpolationError3",
      overwrite = TRUE,
      chunksize = 1e+6
    )    
    
    ## for comparison
    interpolation4 = interpolate(
      simulations = radioactivePlumes_area,
      locations = c(NA, NA, 0, 0, 0, sampleLocations2000, rev(sampleLocations2000)),
      values = 2,
      fun_interpolation = krige0var
    )
    expect_equal(
      names(interpolationError2),
      names(interpolationError4)
    )
    expect_equal(# "cost" (= 'cost_cellStats' as data not in memory -> fun_Rpl cannot be applied)
      interpolationError4[["cost"]],
      cellStats(interpolationError4[["error_locationsplumes"]], "mean")
    )
    
    i = sample.int(nLocations(radioactivePlumes_area), 1)
    j = sample.int(nPlumes(radioactivePlumes_area), 1)  
    expect_equivalent(# "error_locationsplumes"
      delineationError(
        x = c(interpolationError4[["interpolated"]][i,j],
              radioactivePlumes_area@values[[2]][i,j])),
      as.numeric(interpolationError4[["error_locationsplumes"]][i,j,])
    )
    expect_equivalent(# "interpolated"
      interpolationError4[["interpolated"]],
      interpolation4
    )
    
    
})
