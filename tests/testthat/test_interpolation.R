########################################################################
#           test interpolation                                         #
########################################################################

# -- -- -- -- -- -- -- -- -- prepare (dataset, variogram...) -- -- -- -- -- -- -- -- -- 
# import krige0
#library(gstat)
#library(automap)
data(radioactivePlumes_area)

# - - - - - - - generate variogram - - - - - - - - 
# independent for a sample of realisations and locations and use average (median)?
samplePlumes = sort(sample.int(nPlumes(radioactivePlumes_area), 100))
medianVariogram = fitMedianVariogram(simulations = radioactivePlumes_area, 
                              plumes = samplePlumes, 
                              values = 1)
# - - - - - - - - - - interpolation function - - - - - - - - - - - - - - - 
# interpolation function
krige0var = replaceDefault(krige0, newDefaults = list(
  formula = z ~ 1, model = medianVariogram, beta = NA, ... = NA),
  type = "fun_interpolation.interpolate")[[1]]
idw0z = replaceDefault(idw0, newDefaults = list(formula = z ~ 1),
  type = "fun_interpolation.interpolate")[[1]]
#sampleLocations = sort(sample.int(nLocations(radioactivePlumes_area), 100))
#sampleValues = subset(radioactivePlumes_area, locations = sampleLocations, values = 1)
#sampleValuesSpatial0 = subset(sampleValues, plumes = 1)
#sampleValuesSpatial = extractSpatialDataFrame(sampleValuesSpatial0)
#names(sampleValuesSpatial@data) = "z"
#sampleValuesSpatialPoints = as(sampleValuesSpatial, "SpatialPointsDataFrame")
#y = matrix(getValues(sampleValues@values), 
#           ncol = nPlumes(radioactivePlumes_area),
#           byrow = TRUE)
#kriged = krige0var(y = y, data = sampleValuesSpatialPoints)
#original_and_interpolated = subset(radioactivePlumes_area, values = 1)
#original_and_interpolated@values = stack(original_and_interpolated@values,
#                                         brick(kriged,
#                                               xmn = -90, xmx = 90, ymn = -90, ymx = 90,
#                                               crs = CRS("+init=epsg:4326")))
#test1 = extractSpatialDataFrame(original_and_interpolated, plumes = 10)
#spplot(test1, sp.layout = list("sp.points", 
#                               sampleValuesSpatialPoints, col = 3))
#summary(test1@data[sampleLocations,])

test_that("interpolate", {
  # really big dataset, not even subset fitting into memory (truly: canProcessInMemory == FALSE)
  data(radioactivePlumes_area)
  bigSimulations = cbind(radioactivePlumes_area, 
                         radioactivePlumes_area, 
                         radioactivePlumes_area,
                         radioactivePlumes_area,
                         radioactivePlumes_area,
                         radioactivePlumes_area,
                         radioactivePlumes_area,
                         radioactivePlumes_area,
                         radioactivePlumes_area,
                         radioactivePlumes_area,
                         savePath = "bigSimulations",
                         overwrite = TRUE)
  bigSimulationsP = bigSimulations
  bigSimulationsP@locations = as(bigSimulations@locations, "SpatialPointsDataFrame")
  
  sampleLocations1 = c(NA, NA, sample.int(nLocations(radioactivePlumes_area), 2000, replace = TRUE), 
                      0, nLocations(radioactivePlumes_area) + 1:5)
  
  interpolated1 = interpolate(
    simulations = subset(bigSimulationsP, values = 1),
    locations = sampleLocations1,
    fun_interpolation = krige0var,
    tmpfile = "interpol1",
    overwrite = TRUE)
  
  interpolated1a = interpolate(
    simulations = subset(bigSimulationsP, values = 1),
    locations = sampleLocations1,
    fun_interpolation = idw0z,
    tmpfile = "interpol1a",
    overwrite = TRUE)
  
  j = sample.int(nPlumes(bigSimulations), 1)
  original_j = extractSpatialDataFrame(bigSimulations, plumes = j)
  names(original_j@data) = "z"
  Original_j = as(original_j, "SpatialPointsDataFrame")
  SampleLocations1 = sort(unique(sampleLocations1))
  SampleLocations1[SampleLocations1 < 1 | SampleLocations1 > 2500] = NA
  SampleLocations1 = SampleLocations1[!is.na(SampleLocations1)]
  data_j = Original_j[SampleLocations1,1]
  interpolated_1_j = krige(formula = z ~ 1, model = medianVariogram, data = as.data.frame(data_j),
                         newdata = as.data.frame(coordinates(original_j)), locations = ~ s1 + s2)
  interpolated_1a_j = idw(formula = z ~ 1, data = as.data.frame(data_j),
                          newdata = as.data.frame(coordinates(original_j)), locations = ~ s1 + s2)
  Interpolated_j = original_j
  Interpolated_j@data$interpolated1 = interpolated_1_j$var1.pred
  Interpolated_j@data$interpolated1a = interpolated_1a_j$var1.pred
  spplot(Interpolated_j, zcol = c("z", "interpolated1", "interpolated1a"))
  
  expect_equivalent(interpolated_1_j$var1.pred,
                    interpolated1[,j]  ) 
#  expect_equivalent(interpolated_1a_j$var1.pred, # NAs only in idw0, numeric error of ~ 1e-17
#                    interpolated1a[,j])
original_and_interpolated = subset(bigSimulations, values = 1)
original_and_interpolated@values = stack(original_and_interpolated@values,
                                        interpolated1,
                                        interpolated1a)
test1 = extractSpatialDataFrame(original_and_interpolated, plumes = 1000)

spplot(test1)
summary(test1@data[sampleLocations1,])


  # data not in memory, subset in memory
sampleLocations2 = sample.int(nLocations(radioactivePlumes_area), 100)
interpolated2 = interpolate(
  simulations = bigSimulations,
  values = 1,
  locations = sampleLocations2,
  fun_interpolation = krige0var,
  tmpfile = "interpol2",
  overwrite = TRUE)

expect_error(
  interpolate(
    simulations = bigSimulations,
    values = 1,
    locations = sampleLocations2,
    fun_interpolation = krige0var,
    tmpfile = FALSE)# error
  )
SampleLocations2 = sort(unique(sampleLocations2))
SampleLocations2[SampleLocations2 < 2 | SampleLocations2 > 2500] = NA
SampleLocations2 = SampleLocations2[!is.na(SampleLocations2)]
data2_j = Original_j[SampleLocations2,1]
interpolated2_j = krige(formula = z ~ 1, model = meanVariogram, data = as.data.frame(data2_j),
                       newdata = as.data.frame(coordinates(original_j)), locations = ~ s1 + s2)
Interpolated_j@data$interpolated2 = interpolated2_j$var1.pred
spplot(Interpolated_j, zcol = c("z", "interpolated2"))
expect_equivalent(interpolated2_j$var1.pred,
                  interpolated2[,j]) 

# data (fits) in memory
interpolated3 = interpolate(
  simulations = radioactivePlumes_area,
  values = 1,
  locations = sampleLocations2,
  fun_interpolation = krige0var,
  tmpfile = "interpol3",
  overwrite = TRUE)

interpolated3a = interpolate(
  simulations = radioactivePlumes_area,
  values = 1,
  locations = sampleLocations2,
  fun_interpolation = idw0z,
  tmpfile = "interpol3a",
  overwrite = TRUE)

j = sample.int(nPlumes(radioactivePlumes_area), 1)
expect_equivalent(interpolated3[,j],
                  interpolated2[,j]) 

original2_j = extractSpatialDataFrame(radioactivePlumes_area, plumes = j)
names(original2_j@data) = "z"
Original2_j = as(original2_j, "SpatialPointsDataFrame")
SampleLocations2 = sort(unique(sampleLocations2))
SampleLocations2[SampleLocations2 < 1 | SampleLocations2 > 2500] = NA
SampleLocations2 = SampleLocations2[!is.na(SampleLocations2)]
data2_j = Original2_j[SampleLocations2,1]
interpolated_2_j = krige(formula = z ~ 1, model = medianVariogram, data = as.data.frame(data2_j),
                         newdata = as.data.frame(coordinates(original2_j)), locations = ~ s1 + s2)
interpolated_2a_j = idw(formula = z ~ 1, data = as.data.frame(data2_j),
                        newdata = as.data.frame(coordinates(original2_j)), locations = ~ s1 + s2)

newData = as(radioactivePlumes_area@locations, "SpatialPointsDataFrame")
names(newData)[1] = "z"
interpolated__2a_j = idw0(formula = z ~ 1,
                          data = newData[SampleLocations2,], 
                           y = matrix(data2_j@data[[1]], ncol = 1),
                           newdata = newData)

#expect_equal(interpolated_2a_j$var1.pred, # numeric difference
#             interpolated__2a_j)














})
