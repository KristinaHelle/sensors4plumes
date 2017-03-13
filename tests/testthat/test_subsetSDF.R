###############################################################
# test             subsetSDF                                  #
###############################################################
sub_Index = subsetSDF(SIndexDF, locations = c(1,3), data = "c")
sub_Points = subsetSDF(SPointsDF, locations = c(1,3), data = "z")
sub_Pixels = subsetSDF(SPixelsDF, locations = c(1,3), data = "z")
sub_Polygrid = subsetSDF(SPolygridDF, locations = c(1,3), data = "b")
sub_Polygons = subsetSDF(SPolygonsDF, locations = c(1,3), data = "a")
sub_Lines = subsetSDF(SLinesDF, locations = 1:2, data = "a")

test_that("subsetSDF",{
  # class kept
  expect_is(
    sub_Index,
    "SpatialIndexDataFrame"
  )
  expect_is(
    sub_Points,
    "SpatialPointsDataFrame"
  )
  expect_is(
    sub_Pixels,
    "SpatialPixelsDataFrame"
  )
  expect_is(
    sub_Polygrid,
    "SpatialPolygridDataFrame"
  )
  expect_is(
    sub_Polygons,
    "SpatialPolygonsDataFrame"
  )
  expect_is(
    sub_Lines,
    "SpatialLinesDataFrame"
  )
  
  # SpatialIndexDataFrame
  expect_equal(
    sub_Index@data,
    SIndexDF@data[c(1,3), "c", drop = FALSE]
    )
  expect_equal(
    rank(sub_Index@index),
    rank(SIndexDF@index[is.element(SIndexDF@index, c(1,3))])
  )
  ## multiple (ignored) and reordered locations
  grid_index = c(3:7, 12:9)
  locations_index = c(NA,3,4,1,1,3,4,0,1,2)
  sub_Index1 = subsetSDF(SIndexDF, grid = grid_index, locations = locations_index)
  expect_equal(sub_Index1@index, 
               c(2, 2, 2, 2, 1, 1, 1, 1, 2))
  expect_equal(sub_Index1@data, 
               SIndexDF@data[c(3,1),])  
  SIndexDF2 = SpatialIndexDataFrame(
    index = as.integer(c(4,6,1,1,3,2,5,5,5,2,2,3,1,6,2,3,1,4,4,4)),
    data = data.frame(a = 1^{1:6}, b = 0.1 * 1:6, c = c("A", "A", "B", "B", "C", "C"))
    )
  grid_index2 = c(6:1,23,NA,10:20)
  locations_index2 = c(NA,5,5,3,6,3,6,4,1,1,3,4,0,1)
  sub_Index2 = subsetSDF(SIndexDF2, grid = grid_index2, locations = locations_index2)
  sub_Index3 = subsetSDF(SIndexDF2, grid = grid_index2)
  
  expect_equal(
    sub_Index2@index,
    c(3, 2, 4, 4, 1, 1, 4, 2, 1, 4, 3, 3, 3)
  )
  expect_equal(
    sub_Index2@data,
    SIndexDF2@data[c(3,6,4,1),]
  )
  expect_equal(
    sub_Index3@index,
    c(4, 5, 1, 1, 3, 2, 2, 2, 3, 1, 5, 2, 3, 1, 4, 4, 4)
  )
  expect_equal(
    sub_Index3@data,
    SIndexDF2@data[c(1:4, 6),]
  )
  # SpatialPolygridDataFrame
  locations_index2 = c(NA,3,4,1,11,11,1,9,3,9,9,4,0,1,2)
  sub_Polygrid2 = subsetSDF.SpatialPolygridDataFrame(x = SPolygridDF,
                            locations = locations_index2, 
                            coord_x = c(4,12), 
                            coord_y = c(0,8))
  sub_Polygrid3 = subsetSDF.SpatialPolygridDataFrame(x = SPolygridDF,
                            locations = locations_index2, 
                            grid_i = c(3,6,5,4), 
                            grid_j = c(3,5,6,4))                           
  sub_Polygrid4 = subsetSDF.SpatialPolygridDataFrame(x = SPolygridDF,
                            locations = locations_index2, 
                            grid_ij = matrix(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                               FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                               FALSE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE,
                                               FALSE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE,
                                               FALSE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE,
                                               FALSE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE),
                                             nrow = 6, byrow = TRUE))
  expect_equal(
    sub_Polygrid2@index,
    c(3, 6, 5, 5, 2, 1, 5, 5, 4, 4, NA, NA, 4, 4, NA, NA)
  )
  expect_equal(
    sub_Polygrid2@data,
    SPolygridDF@data[c(3,4,1,11,9,2),]
  )
  expect_equal(
    sub_Polygrid2,
    sub_Polygrid3
  )
  expect_equal(
    sub_Polygrid2,
    sub_Polygrid4
  )
  # from sp
  expect_equal(
    sub_Points,
    SPointsDF[c(1,3), "z"]
  )
  expect_equal(
    sub_Pixels,
    SPixelsDF[c(1,3), "z"]
  )
  expect_equal(
    sub_Pixels,
    SPixelsDF[c(1,3), "z"]
  )
  expect_equal(
    sub_Polygons,
    SPolygonsDF[c(1,3), "a"]
  )
  expect_equal(
    sub_Lines,
    SLinesDF[1:2, "a"]
  )
})
