#############################################
# test         summaryLocations             #
#############################################

test_that("summaryLocations", 
{
  data(radioactivePlumes_area)
  i = sample(nLocations(radioactivePlumes_area), 1)
  # basic functionality with pedefined fun
  summaryLocations_1 = summaryLocations(radioactivePlumes_area, fun = sum)
  expect_equivalent(
    summaryLocations_1[[2]][i],
    sum(getValues(radioactivePlumes_area@values, i, 1)[,1])
  )
  # user-defined fun
  fun1 = function(x, ...){
    xExceed = x > 1e-7
    out = sum(xExceed)
    return(out)
  }
  summaryLocations_2 = summaryLocations(radioactivePlumes_area, fun = fun1, values = 1)
  expect_equivalent(
    summaryLocations_2[[2]][i],
    fun1(getValues(radioactivePlumes_area@values, i, 1)[,1])
  )
  # values
  summaryLocations_3 = summaryLocations(radioactivePlumes_area, 
                                        fun = max, values = "maxdose")  
  expect_equivalent(
    summaryLocations_3[[2]][i],
    max(getValues(radioactivePlumes_area@values, i, 1)[,"maxdose"])
  )
  # na.rm
  summaryLocations_4 = summaryLocations(radioactivePlumes_area, 
                                        fun = max, values = 3, na.rm = TRUE) 
  expect_true(
    all(!is.na(summaryLocations_4[[2]]))
  )
  # (this test may fail for other data if there is no NA in "time")
  summaryLocations_5 = summaryLocations(radioactivePlumes_area, 
                                        fun = max, values = 3, na.rm = FALSE) 
  expect_true(
    sum(is.na(summaryLocations_5[[2]])) > 0
  )
  
  # indices
  plumes = sample(nPlumes(radioactivePlumes_area), 100)
  locations = sample(nLocations(radioactivePlumes_area), 100)
  # locations
  # locations out of bounds -> NA-values
  summaryLocations_6 = summaryLocations(radioactivePlumes_area, fun = sum, 
                                        locations = c(locations, -2, 0, 3000, NA), values = 1, na.rm = FALSE)
  expect_true(
    all(is.na(summaryLocations_6[[2]][101:104]))
  )    
  expect_equivalent(
    as.vector(summaryLocations_6[[2]][1:100]),
    as.vector(summaryLocations_1[[2]][locations])
  )
  # plumes
  summaryLocations_7 = summaryLocations(radioactivePlumes_area, 
                                        fun = sum, plumes = plumes, values = 1)
  summaryLocations_8 = summaryLocations(subset(radioactivePlumes_area, plumes = plumes), 
                                        fun = sum, values = 1)
  expect_equal(
    summaryLocations_7,
    summaryLocations_8
  )
  # plumes out of bounds -> NA at these indices
  expect_warning(summaryLocations_9 <- summaryLocations(radioactivePlumes_area, fun = sum, values = 1,
                                                      plumes = c(plumes, -2, 0, 2324, NA, plumes), na.rm = FALSE))
  expect_true(
    all(is.na(summaryLocations_9[[2]]))  
  )
  expect_warning(summaryLocations_10 <- summaryLocations(radioactivePlumes_area, fun = sum, values = 1,
                                                        plumes = c(plumes, -2), na.rm = TRUE))
  expect_warning(summaryLocations_11 <- summaryLocations(radioactivePlumes_area, fun = sum, values = 1,
                                                         plumes = plumes, na.rm = TRUE))
  expect_equivalent(
    summaryLocations_10,
    summaryLocations_11
  )  

  # weight
  # weight (numeric)
  summaryLocations_12 = summaryLocations(radioactivePlumes_area, fun = sum, values = 1, 
                                         weight = 1:nLocations(radioactivePlumes_area))
  expect_equal(
    summaryLocations_12[[1]],
    mean(summaryLocations_1[[2]] * 1:nLocations(radioactivePlumes_area))
  )

  radioactivePlumes_area@locations@data$dist2origin = 
    spDistsN1(coordinates(radioactivePlumes_area@locations), coordinates(radioactivePlumes_area@locations)[1,])
  
  # weight (character)
  summaryLocations_13 = summaryLocations(radioactivePlumes_area, fun = sum, values = 1, 
                                       weight ="dist2origin")
  expect_equal(
    summaryLocations_13[[1]],
    mean(summaryLocations_13[[2]] * radioactivePlumes_area@locations@data$dist2origin) 
  )
  
  # weight (character), locations
  weightedMeanNArm = replaceDefault(weightedMean, newDefaults = list(na.rm = TRUE), 
                                    type = "summaryFun.summaryPlumes")
  expect_warning(summaryLocations_14 <- summaryLocations(radioactivePlumes_area, fun = sum, kinds = 1, 
                                                         locations = c(NA, locations),
                                         weight ="dist2origin", summaryFun = weightedMeanNArm[[1]]))
  expect_equal(
    summaryLocations_14[[1]],
    mean(summaryLocations_14[[2]] * c(NA, radioactivePlumes_area@locations@data$dist2origin[locations]), na.rm = TRUE) 
  )
})

