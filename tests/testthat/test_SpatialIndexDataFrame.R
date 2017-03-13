#######################################################################
# test        SpatialIndexDataFrame                                   #
#######################################################################
# SpatialIndexDataFrame()
SIndexDF1 = SpatialIndexDataFrame(index1, dataFrame1) 
SIndexDF2 = SpatialIndexDataFrame(as.integer(1), dataFrame1[1,])
test_that("SpatialIndexDataFrame()", {
  expect_error(# wrong order
    SpatialIndexDataFrame(dataFrame1, index1) 
  )  
  expect_warning(# numbers in index do not fit rows of data
    SIndexDF3 <-SpatialIndexDataFrame(index = as.integer(c(5,5,NA,3,4,6,6,3,-7,18,4,3)), dataFrame1)
  )
  expect_equal(
    SIndexDF3@data,
    dataFrame1[3:6,]
  )
  expect_equal(
    SIndexDF3@index,
    as.integer(c(3,3,NA,1,2,4,4,1,NA,NA,2,1))
  )
})

test_that("coerce:data.frame -> SpatialIndexDataFrame", {
  SIndexDF3 = as(dataFrame3, "SpatialIndexDataFrame")
  expect_equal(
    SIndexDF3@data,
    dataFrame3
  )
  expect_equal(
    SIndexDF3@index,
    1:nrow(dataFrame3)
  )
})

test_that("coordinates.SpatialIndexDataFrame", {
  expect_warning(
    coordinatesB1 <- coordinates(SIndexDF)  
  )
  expect_null(
    coordinatesB1
  )
})
test_that("bbox.SpatialIndexDataFrame", {
  expect_warning(
    bboxB1 <- bbox(SIndexDF)     
  )
  expect_equal(
    bboxB1,
    matrix(0, nrow = 2, ncol = 2)
  )
})
test_that("proj4string.SpatialIndexDataFrame", {
  expect_warning(
    proj4stringB1 <- proj4string(SIndexDF)
  )
  expect_equal(
    proj4stringB1,
    as.character(NA)
  )
})
test_that("spplot.SpatialIndexDataFrame", {
  expect_warning(
    spplot(SIndexDF)
  )  
  expect_error(
    spplot(SIndexDF2)
  )
})

test_that("cbind.SpatialIndexDataFrame", {
  dataFrameA1 = data.frame(f = c(10, 20, 30), 
                                  g = c(1,2,4))
  SIndexDF4 = cbind(SIndexDF, 
                    SpatialIndexDataFrame(index = SIndexDF@index,
                                          data = dataFrameA1)
                    )
  expect_equal(
    SIndexDF4@data,
    cbind(SIndexDF@data, dataFrameA1)
  )
  expect_equal(
    SIndexDF4@index,
    SIndexDF@index
  )
  expect_error(
    cbind(SIndexDF,SIndexDF3)
  )
})
test_that("rbind.SpatialIndexDataFrame", {
  SIndexDF5 = rbind(SIndexDF1, SIndexDF2)
  expect_equal(
    SIndexDF5@index,
    c(SIndexDF1@index, 13)
  )  
  expect_equal(
    dename(SIndexDF5@data, kind = "rownames"),
    dename(dataFrame1[c(1:12,1),], kind = "rownames")
  ) 
  expect_error(
    rbind(SIndexDF, SIndexDF1)  
  )
})

# subsetSDF => subsetSDF
# areaSDF => areaSDF