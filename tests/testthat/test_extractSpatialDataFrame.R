#########################################################
# test extractSpatialDataFrame                          #
#########################################################

# examples of all SDF classes (needs objects of test_Simulations.R, test_pbind.R)
eSIndexDF1 = extractSpatialDataFrame(SimulationsIndex1)
eSPointsDF1 = extractSpatialDataFrame(SimulationsPoints1)
eSPixelsDF1 = extractSpatialDataFrame(SimulationsPixels1)
eSPolygridDF1 = extractSpatialDataFrame(SimulationsPolygrid1)
eSPolygonsDF1 = extractSpatialDataFrame(SimulationsPolygons1)
eSLinesDF1 = extractSpatialDataFrame(SimulationsLines1)

# small example
data(SimulationsSmall)
i = sample.int(nPlumes(SimulationsSmall), 4)
j = sample.int(nKinds(SimulationsSmall), 2)
SimulationsSmall_ij = extractSpatialDataFrame(
  SimulationsSmall,
  plumes = i, kinds = j)

# big example
data(radioactivePlumesArea)
RadioactivePlumes_area_725_maxdose = extractSpatialDataFrame(
  RadioactivePlumes_area, 
  plumes = c(725, 15), kinds = c(3,2))

# example too big to load all 
expect_warning(
  RadioactivePlumes_area2 <- extractSpatialDataFrame(
    RadioactivePlumes_area))
expect_warning(
  RadioactivePlumes_area3 <- extractSpatialDataFrame(
    RadioactivePlumes_area, plumes = 800:1, kinds = 3:1))
                                                             
test_that("extractSpatialDataFrame", {
  # correct classes
  expect_is(eSIndexDF1, "SpatialIndexDataFrame")
  expect_is(eSPointsDF1, "SpatialPointsDataFrame")
  expect_is(eSPixelsDF1, "SpatialPixelsDataFrame")
  expect_is(eSPolygridDF1, "SpatialPolygridDataFrame")
  expect_is(eSPolygonsDF1, "SpatialPolygonsDataFrame")
  expect_is(eSLinesDF1, "SpatialLinesDataFrame")
  
  expect_equal(
    proj4string(eSPointsDF1),
    proj4string(SimulationsPoints1@locations))
  
  # correct values
  k = sample.int(length(SimulationsSmall_ij), 1)
  expect_equal(as.numeric(SimulationsSmall_ij@data[k,]),
               as.numeric(SimulationsSmall@values[[j]][k,i]))
  k = sample.int(length(RadioactivePlumes_area_725_maxdose), 1)
  expect_equal(as.numeric(RadioactivePlumes_area_725_maxdose@data[k,]),
               as.numeric(getValues(subset(RadioactivePlumes_area@values, c(3,2)), row = k, nrows = 1)[c(725, 15),]))

  expect_equal(as.numeric(RadioactivePlumes_area2@data[k,]),
               as.numeric(getValues(subset(RadioactivePlumes_area@values, 1), 
                                    row = k, nrows = 1)[1:1600]))
  expect_equal(as.numeric(RadioactivePlumes_area3@data[k,]),
               as.numeric(getValues(subset(RadioactivePlumes_area@values, c(3,2)), 
                                    row = k, nrows = 1)[800:1,]))
  
  # errors and warnings (too many, invalid, or no valid indices)  
 expect_warning(
   SimulationsSmall_32 <- extractSpatialDataFrame(SimulationsSmall, plume = c(6,3,1), kind = 2))
 expect_warning(
   expect_error( 
     extractSpatialDataFrame(SimulationsSmall, plumes = 2, kinds = 3)))
 expect_warning(
   expect_error( 
     extractSpatialDataFrame(SimulationsSmall, plumes = c(-1, 0), kinds = 1)))
})


  
