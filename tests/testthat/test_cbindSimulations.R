#########################################################################
# test                        cbind (Simulations)                       #
#########################################################################
# needs all data from test_Simulations.R

# locations: all kind of SDF
SimulationsIndex2 = cbind(SimulationsIndex1, SimulationsIndex1)
test_that("cbind SIndexDF", {
  nP = nPlumes(SimulationsIndex1)
  expect_is(SimulationsIndex2, "Simulations")
  expect_equal(nPlumes(SimulationsIndex2), 2 * nP)
  expect_equal(nKinds(SimulationsIndex2), nKinds(SimulationsIndex1))
  expect_equal(SimulationsIndex2@locations, SimulationsIndex1@locations)
  expect_equal(SimulationsIndex2@plumes, 
               rbind(SimulationsIndex1@plumes, SimulationsIndex1@plumes))
  expect_equal(SimulationsIndex2@values[,][nP + 1:nP,], 
               SimulationsIndex1@values[,][1:nP,])
})

SimulationsPoints2 = cbind(SimulationsPoints1, SimulationsPoints1)
SimulationsPixels2 = cbind(SimulationsPixels1, SimulationsPixels1)
SimulationsPolygrid2 = cbind(SimulationsPolygrid1, SimulationsPolygrid1)
SimulationsPolygons2 = cbind(SimulationsPolygons1, SimulationsPolygons1)
SimulationsLines2 = cbind(SimulationsLines1, SimulationsLines1)

SimulationsPoints3 = cbind(SimulationsPoints1, SimulationsPoints1, SimulationsPoints1)

test_that("cbind", {
  expect_equal(nPlumes(SimulationsPoints2), 2 * nPlumes(SimulationsPoints1))
  expect_equal(nPlumes(SimulationsPoints3), 3 * nPlumes(SimulationsPoints1))
  expect_equal(nPlumes(SimulationsPixels2), 2 * nPlumes(SimulationsPixels1))
  expect_equal(nPlumes(SimulationsPolygrid2), 2 * nPlumes(SimulationsPolygrid1))
  expect_equal(nPlumes(SimulationsPolygons2), 2 * nPlumes(SimulationsPolygons1))
  expect_equal(nPlumes(SimulationsLines2), 2 * nPlumes(SimulationsLines1))
})

# big example
if (FALSE){
  data(radioactivePlumesArea) # adapt to final name of file
  RadioactivePlumes_area2 = cbind(RadioactivePlumes_area, RadioactivePlumes_area,
                                  savePath = "/home/kristina/Desktop/s4p_article/test/cbind1_")
  test_that("cbind big SPixelsDF", {
    nP = nPlumes(RadioactivePlumes_area)
    expect_is(RadioactivePlumes_area2, "Simulations")
    expect_equal(nPlumes(RadioactivePlumes_area2), 2 * nP)
    expect_equal(nKinds(RadioactivePlumes_area2), nKinds(RadioactivePlumes_area))
    expect_equal(RadioactivePlumes_area2@locations, RadioactivePlumes_area@locations)
    expect_equal(RadioactivePlumes_area2@plumes, 
                 rbind(RadioactivePlumes_area@plumes, RadioactivePlumes_area@plumes))
    expect_equal(getValues(RadioactivePlumes_area2@values, 1, 1)[nP + 1:nP,], 
                 getValues(RadioactivePlumes_area@values, 1, 1)[1:nP,])
  })  
}


# non-fitting examples: 
## locations differ by points
SimulationsPoints1A = Simulations(locations = subsetSDF(SPointsDF, locations = 3:4),
                                 plumes = plumes1,
                                 values = values1)
## locations differ by projection
SimulationsPoints1a = SimulationsPoints1
SimulationsPoints1b = SimulationsPoints1
proj4string(SimulationsPoints1a@locations) = CRS("+init=epsg:4267")
proj4string(SimulationsPoints1b@locations) = CRS("+init=epsg:26978")
## locations differ by attributes
SimulationsPoints1B = SimulationsPoints1
SimulationsPoints1B@locations@data$z = 3:4

## plumes differ
### different attribute name
SimulationsPixels1A = SimulationsPixels1
names(SimulationsPixels1A@plumes)[1] = "Source"
### different attribute type
SimulationsPixels1B = SimulationsPixels1
SimulationsPixels1B@plumes$source = as.character(SimulationsPixels1@plumes$source) 
SimulationsPixels2B = cbind(SimulationsPixels1, SimulationsPixels1B) # works, more complex type is used for both

## value layers differ
### differ by value, should not be a problem
values1A = stack(raster(x = replicate(n = 5, expr = rlnorm(2, sdlog = 2)), 
                       xmn = -90, xmx = 90, ymn = -90, ymx = 90, 
                       crs = "+init=epsg:4326" ), 
                raster(x = replicate(n = 5, expr = rnorm(2, m = 10, sd = 3)), 
                       xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                       crs = "+init=epsg:4326" ))
SimulationsPolygrid1A = Simulations(locations = subsetSDF(SPolygridDF, locations = 1:2),
                                    plumes = plumes1,
                                    values = values1A)
SimulationsPolygrid3 = cbind(SimulationsPolygrid1, SimulationsPolygrid1A)
### differ by dataType
raster1 = raster(x = replicate(n = 5, expr = rlnorm(2, sdlog = 2)), 
                 xmn = -90, xmx = 90, ymn = -90, ymx = 90, 
                 crs = "+init=epsg:4326" )
raster2 = raster(x = replicate(n = 5, expr = rbinom(2, 1, 0.5)), 
                 xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                 crs = "+init=epsg:4326" )
dataType(raster2) = "LOG1S"
values2 = stack(raster1, raster2)
SimulationsPolygrid1B = Simulations(locations = subsetSDF(SPolygridDF, locations = 1:2),
                                    plumes = plumes1,
                                    values = values2)
SimulationsPolygrid4 = cbind(SimulationsPolygrid1, SimulationsPolygrid1B) # works, more complex type is used for both

### differ by name
SimulationsPolygrid1C = SimulationsPolygrid1
names(SimulationsPolygrid1C@values)[1] = "Layer.1"

test_that("cbind errors", {
  # locations differ
  expect_error(
    cbind(SimulationsPoints1, SimulationsPoints1A))
  expect_error(
    cbind(SimulationsPoints1a, SimulationsPoints1b)) 
  expect_error(
    cbind(SimulationsPoints1, SimulationsPoints1B))
  # plumes differ
  expect_error(
    cbind(SimulationsPixels1, SimulationsPixels1A))
  # values differ (by name)
  expect_error(
    cbind(SimulationsPolygrid1, SimulationsPolygrid1, SimulationsPolygrid1C))
    
  # no error (more complex data type is used for combination)
  expect_is(SimulationsPixels2B@plumes$source, "factor")
  expect_equal(dataType(SimulationsPolygrid4@values)[2], "FLT4S")
})

## differences between locations (locations, proj4string..)
## differences between values (number of layers), name of layers, datatype of layers (?)
## differences between plumes: number, name, type

