###################################################################
# test SpatialPolygridDataFrame()                                 #
###################################################################
SPolygridDFA1 = SpatialPolygridDataFrame(
  grid = grid1,
  data = dataFrame1,
  index = index1)
coordinatesA1 = coordinates(SPolygridDFA1)

SPolygridDFA2 = SpatialPolygridDataFrame(
  data = dataFrame2,
  index = index1,
  grid = grid1,
  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

SPolygridDFA3 = SpatialPolygridDataFrame(
  data = dataFrame2,
  index = index1,
  grid = spatialGrid1)

test_that("SpatialPolygridDataFrame()", {
  expect_is(SPolygridDFA1, "SpatialPolygridDataFrame")
  expect_equal(SPolygridDFA2, SPolygridDFA3)
  
  expect_error( # length of index does not fit grid
    SpatialPolygridDataFrame(
      grid = grid1,
      data = dataFrame1,
      index = as.integer(c( 6, 7, 7, 8, 8,
                            6, 6, 7, 7, 8, 8,
                            5, 5, 1, 2, 9, 9,
                            5, 5, 4, 3, 9, 9,
                            12,12,11,11,10,10,
                            12,12,11,11,10,10)))
  )
  
  expect_warning( # numbers in index do not fit length of data
    SPolygridDFF4 <- SpatialPolygridDataFrame(
      grid = grid1,
      data = dataFrame1,
      index = as.integer(c( 6, 6, 7, 7, 8, 8,
                            6, 6, 7, 7, 8, 8,
                            5, 5, 1, 2, 9, 9,
                            5, 5, 4, 3, 9, 9,
                            11,11,11,11,10,10,
                            11,11,11,11,10,10)))
  )
  expect_equal(
    SPolygridDFF4@data,
    dataFrame1[1:11,]
  )
  
  expect_warning( # numbers in index not complete sequence from 1
    SPolygridDFF5 <- SpatialPolygridDataFrame(
      grid = grid1,
      data = dataFrame1,
      index = as.integer(c( 6, 6, 7, 7, -8, -8,
                            6, 6, 7, 7, -8, -8,
                            5, 5, 13, 2, 9, 9,
                            5, 5, 4, 3, 9, 9,
                            -12,-12,11,11,10,10,
                            -12,-12,11,11,10,10)))
  )
  expect_equal(
    SPolygridDFF5@data,
    dataFrame1[c(2:7, 9:11),]
  )
  expect_equal(
    SPolygridDFF5@index,
    as.integer(c( 5, 5, 6, 6,NA,NA,
                  5, 5, 6, 6,NA,NA,
                  4, 4,NA, 1, 7, 7,
                  4, 4, 3, 2, 7, 7,
                  NA,NA, 9, 9, 8, 8,
                  NA,NA, 9, 9, 8, 8))
  )  
  
  
  expect_error( # parameters in wrong order data <-> index
    SpatialPolygridDataFrame(
      index1,
      dataFrame1,
      grid1)
  )
  
  expect_error( # parameters in wrong order grid <-> index
    SpatialPolygridDataFrame(
      dataFrame1,
      grid1,
      index1)
  )
  
  expect_error( # parameters in wrong order grid <-> data
    SpatialPolygridDataFrame(
      grid1,
      index1,
      dataFrame1)
  )
  
  expect_error(
    SpatialPolygridDataFrame( # invalid proj4string argument: CRS missing
      dataFrame1,
      index1,
      grid1,
      "+proj=longlat")
  )
  
  expect_error(
    SpatialPolygridDataFrame( # invalid proj4string argument in CRS
      dataFrame1,
      index1,
      grid1,
      CRS("+ proj = longlat"))
  )    
})  
