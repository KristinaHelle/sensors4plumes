##############################################################
# test  SpatialDataFrame    coordinates  proj4string  bbox   #
##############################################################
coordinates3 = coordinates(SPolygridDFA1)
# --------------------- coordinates ------------------------ #
# when transforming regular points into SPolygridDF, the resulting coordinates should be the same
test_that("coordinates.SpatialDataFrame", {
  i = sample.int(12,1)
  expect_equal(
    coordinates(SPolygridDFA1)[i,], 
    apply(spatialPoints1@coords[SPolygridDFA1@index == i,,drop = FALSE], mean, MARGIN = 2))
})

# --------------------- proj4string ------------------------ #
test_that("proj4string.SpatialDataFrame", {
  expect_equal(proj4string(SPolygridDFB2), 
               "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  )
  expect_error( # proj4string can only retrieve value, not set it
    proj4string(SPolygridDFB2) <- CRS("+proj=longlat +datum=WGS84")
  )    
})

# --------------------- bbox -------------------------------- #
test_that("bbox.SpatialDataFrame", {
  expect_equal(
    dename(bbox(SPolygridDFA1), kind = "dimnames"),
    matrix(c(0, 0, 12, 12), nrow = 2)   
  )  
})
