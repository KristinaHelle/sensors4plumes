#############################################
# test         summaryPlumes                # 
#############################################

### user-defined fun: it cannot be tested beforehand if function is truly associative and only returns one value
test_that("summaryPlumes",{
  
# big example (not in memory)
data(radioactivePlumes_area)
i = sample.int(nPlumes(radioactivePlumes_area), 1) 
# proper use of indices
plumes = sample.int(nPlumes(radioactivePlumes_area), 100)
locations = sample.int(nLocations(radioactivePlumes_area), 100)

## basic functionality with predefined fun
summaryPlumes_1 = summaryPlumes(radioactivePlumes_area, fun = sum, kinds = 2)
summaryPlumes_2 = summaryPlumes(radioactivePlumes_area, fun = prod, kinds = 2)
summaryPlumes_3 = summaryPlumes(radioactivePlumes_area, fun = min, kinds = 2)
summaryPlumes_4 = summaryPlumes(radioactivePlumes_area, fun = max, kinds = 2)
expect_equivalent(
  summaryPlumes_1[[2]][i],
  sum(radioactivePlumes_area@values[,i][,2])
) 
expect_equivalent(
  summaryPlumes_3[[2]][i],
  min(radioactivePlumes_area@values[,i][,2])
)
expect_equivalent(
  summaryPlumes_4[[2]][i],
  max(radioactivePlumes_area@values[,i][,2])
)
## effect of na.rm
summaryPlumes_5 = summaryPlumes(radioactivePlumes_area, fun = sum, kinds = 3)
summaryPlumes_6 = summaryPlumes(radioactivePlumes_area, fun = sum, kinds = 3, na.rm = TRUE)
expect_true(
  sum(is.na(summaryPlumes_5[[2]])) > 0
) 
expect_true(
  sum(is.na(summaryPlumes_6[[2]])) == 0
) 
## user-defined fun (correct)
summaryPlumes_7 = summaryPlumes(radioactivePlumes_area, 
                                fun = function(x, ...){sum(x, ...)}, kinds = 2)
expect_equal(
  summaryPlumes_1,
  summaryPlumes_7
)
## user-defined fun (incorrect: missing ... or na.rm)
expect_error(
  summaryPlumes(radioactivePlumes_area, 
                fun = function(x){x}, initial = 0, kinds = 2)
)  

# include subsetting
# nonexistent locations are counted NA
expect_warning(
  summaryPlumes_8 <- summaryPlumes(radioactivePlumes_area, fun = sum, 
                        locations = c(locations, -2, 0, 3000), kinds = 2)
)
expect_true(
  all(is.na(summaryPlumes_8[[2]]))
)
expect_warning(
  summaryPlumes_9 <- summaryPlumes(radioactivePlumes_area, fun = sum, na.rm = TRUE,
                                   locations = c(1:2500, -2, 0, 3000), kinds = 2)
)
expect_equal(
  summaryPlumes_1[[1]],
  summaryPlumes_9[[1]]
)
expect_equal(
  summaryPlumes_1[[2]],
  summaryPlumes_9[[2]]
)
# multiple locations are deleted
expect_warning(
  summaryPlumes_10 <- summaryPlumes(radioactivePlumes_area, fun = sum, na.rm = TRUE,
                                   locations = c(locations, -2, 0, 3000, NA, locations), kinds = 2)
)
expect_equivalent(
  summaryPlumes_10[[2]][i],
  sum(radioactivePlumes_area@values[locations,i][,2])
)

# plumes out of bounds -> NA at these indices
# multiple plumes taken into account
expect_warning(
  summaryPlumes_11 <- summaryPlumes(radioactivePlumes_area, fun = sum, 
                                   plumes = c(plumes, -2, 0, 2324, NA, plumes), 
                                   na.rm = TRUE, kinds = 2)
)
expect_equivalent(
  summaryPlumes_11[[2]][10],
  sum(radioactivePlumes_area@values[,plumes[10]][,2])
)
expect_equivalent(
  summaryPlumes_11[[2]][114],
  sum(radioactivePlumes_area@values[,plumes[10]][,2])
)  
expect_true(
  all(is.na(summaryPlumes_11[[2]][101:104]))
)  
expect_true(
  is.na(summaryPlumes_11[[1]]) # because default in weightedMean
)
weightedMeanNArm = replaceDefault(weightedMean, newDefaults = list(na.rm = TRUE), 
                                  type = "summaryFun.summaryPlumes")[[1]]
expect_warning(
  summaryPlumes_12 <- summaryPlumes(radioactivePlumes_area, fun = sum, 
                                    summaryFun = weightedMeanNArm,
                                    plumes = c(plumes, -2, 0, 2324, NA, plumes), 
                                    na.rm = TRUE, kinds = 2)
)
expect_true(  
  !is.na(summaryPlumes_12[[1]]) 
)  
expect_equal(
  summaryPlumes_11[[2]],
  summaryPlumes_12[[2]]
)
  
# weights
# weight (numeric)
summaryPlumes_13 = summaryPlumes(radioactivePlumes_area, fun = sum, 
                                 weight = 1:nPlumes(radioactivePlumes_area), kinds = 2)
expect_equal(
  summaryPlumes_13[[1]],
  mean(summaryPlumes_1[[2]] * 1:nPlumes(radioactivePlumes_area))
)

# weight (character)
radioactivePlumes_area@plumes$totalDose = 
  summaryPlumes(radioactivePlumes_area, kinds = "finaldose", fun = sum)[[2]]
summaryPlumes_14 = summaryPlumes(radioactivePlumes_area, fun = sum, 
                                 weight = "totalDose", kinds = 2)
summaryPlumes_15 = summaryPlumes(radioactivePlumes_area, fun = sum, 
                                 weight = radioactivePlumes_area@plumes$totalDose, kinds = 2)
expect_equal(
  summaryPlumes_14,
  summaryPlumes_15
)  
# weight (character), plumes
expect_warning(
  summaryPlumes_16 <- summaryPlumes(radioactivePlumes_area, fun = sum, weight = "totalDose", 
                                   kinds = 2, plumes = c(plumes, NA, plumes))  
)
expect_warning(
  summaryPlumes_17 <- summaryPlumes(radioactivePlumes_area, fun = sum, 
                                    weight = radioactivePlumes_area@plumes$totalDose[c(plumes, NA, plumes)], 
                                    kinds = 2, plumes = c(plumes, NA, plumes))  
)
expect_equal(
  summaryPlumes_16,
  summaryPlumes_17
)  
}
)
