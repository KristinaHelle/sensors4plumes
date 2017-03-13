##############################################################
# test SpatialPolygridDataFrame coerce                       #
##############################################################
SPolygridDFB1 = as(spatialGridDataFrame1, "SpatialPolygridDataFrame")
spatialGridDataFrameB1 = as(SPolygridDFB1, "SpatialGridDataFrame")
SPolygridDFB2 = as(spatialPointsDataFrame1, "SpatialPolygridDataFrame")
SPolygridDFB3 = as(spatialPointsDataFrame2, "SpatialPolygridDataFrame")
spatialGridDataFrameB2 = as(SPolygridDFB3, "SpatialGridDataFrame")


test_that("coerce: SpatialGridDataFrame -> SpatialPolygridDataFrame", {
  expect_equal(
    SPolygridDFB1@grid, SPolygridDFA2@grid)
  expect_equivalent(
    SPolygridDFB1@data, 
    SPolygridDFA1@data[index1,])
  expect_equivalent(
    SPolygridDFB1@data, 
    SPolygridDFA1@data[index1,])
})

test_that("coerce: SpatialPolygridDataFrame -> SpatialGridDataFrame", {
  expect_equivalent(
    spatialGridDataFrame1@grid@cellcentre.offset, 
    spatialGridDataFrameB1@grid@cellcentre.offset)
  expect_equivalent(
    spatialGridDataFrame1@grid@cells.dim, 
    spatialGridDataFrameB1@grid@cells.dim)
  expect_equivalent(
    spatialGridDataFrame1@grid@cellsize, 
    spatialGridDataFrameB1@grid@cellsize)
  expect_equivalent(
    spatialGridDataFrame1@data, 
    spatialGridDataFrameB1@data)  
})

test_that("coerce: SpatialPointsDataFrame {regular} -> SpatialGridDataFrame", {
  expect_equal(
    SPolygridDFB1, 
    SPolygridDFB2)
})

test_that("coerce: SpatialPointsDataFrame {irregular} -> SpatialGridDataFrame", {
  # bounding box covers all points
  tol = mean(diff(t(bbox(spatialPointsDataFrame2))))/1000
  expect_more_than(SPolygridDFB3@bbox[1,2],
                   spatialPointsDataFrame2@bbox[1,2] - 0.5 * tol)
  expect_more_than(SPolygridDFB3@bbox[2,2],
                   spatialPointsDataFrame2@bbox[2,2] - 0.5 * tol)
  expect_less_than(SPolygridDFB3@bbox[1,1],
                   spatialPointsDataFrame2@bbox[1,1] + 0.5 * tol)
  expect_less_than(SPolygridDFB3@bbox[2,1],
                   spatialPointsDataFrame2@bbox[2,1] + 0.5 * tol)  
  # arbitrary cell has correct value 
  i = sample.int(SPolygridDFB3@grid@cells.dim[1], 1)
  j = sample.int(SPolygridDFB3@grid@cells.dim[2], 1)
  neighbour_ij = which.min(spDistsN1(coordinates2,
                                     spatialGridDataFrameB2[i,j]@coords))
  expect_equivalent(
    spatialPointsDataFrame2@data[neighbour_ij,], 
    spatialGridDataFrameB2[i,j]@data
  )
#  Col = rep(3, length(spatialPointsDataFrame2))
#  Col[neighbour_ij] = "white"
#  spplot(spatialGridDataFrameB2, zcol = "b",
#         sp.layout = list("sp.points", rbind(spatialGridDataFrameB2[i,j]@coords,
#                                             coordinates(spatialPointsDataFrame2)
#                                             ),
#                          col = c(2, Col)))

  # all data occur
  expect_equivalent(
    unique(spatialPointsDataFrame2@data),
    unique(SPolygridDFB3@data)
  )
})

test_that("SpatialPolygridDataFrame -> SpatialPointsDataFrame",{
  data(SPolygridDF)
  SPointsDF1 = as(SPolygridDF, "SpatialPointsDataFrame")
  i = sample.int(length(SPointsDF1), 1)
  SPolygridDF_i = subsetSDF(SPolygridDF, locations = i)
  SPointsDF_i = SPointsDF1[i,]
  expect_equal(
    coordinates(SPolygridDF_i),
    coordinates(SPointsDF_i)
  )
  expect_equal(
    SPolygridDF_i@data,
    SPointsDF_i@data
  )  
})

