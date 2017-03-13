#####################################################################
# test SpatialPolygridDataFrame  cbind                              #
#####################################################################
SPolygridDFE1 = cbind(SPolygridDFA2, subsetSDF(SPolygridDFA3, data = "a"))

test_that("cbind.SpatialPolygridDataFrame", {
  # data is combined
  expect_equal(
    dename(SPolygridDFE1@data),
    dename(dataFrame2[,c("a", "b", "a")])
  )
  # spatial reference is kept
  expect_equal(
    SPolygridDFE1@index, SPolygridDFA2@index
  )
  expect_equal(
    SPolygridDFE1@grid, SPolygridDFA3@grid
  )
  
  # grids do not agree (proj4string differs)
  expect_error(
    cbind(SPolygridDFA1, SPolygridDFA2)
  )
  # indices do not agree
  expect_error(
    cbind(SPolygridDFA2, SPolygridDFB1)
  )  
})


