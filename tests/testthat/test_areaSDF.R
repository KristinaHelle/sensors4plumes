###############################################################
# test             areaSDF                                  #
###############################################################
expect_warning(
  area_Index <- areaSDF(SIndexDF)  
)
expect_warning(
  area_Points <- areaSDF(SPointsDF)
)
area_Pixels = areaSDF(SPixelsDF)
area_Polygrid = areaSDF(SPolygridDF)
area_Polygons = areaSDF(SPolygonsDF)
area_Lines = areaSDF(SLinesDF)

test_that("areaSDF",{
  expect_equal(
    area_Index,
    rep(0, nrow(SIndexDF@data))
  )
  expect_equal(
    area_Points,
    rep(0, nrow(SPointsDF@data))
  )
  expect_equal(
    area_Pixels,
    rep(prod(SPixelsDF@grid@cellsize), nrow(SPixelsDF@data))
  )  
  expect_equal(
    area_Polygrid,
    prod(SPolygridDF@grid@cellsize) * as.data.frame(table(SPolygridDF@index))$Freq
  )   
  expect_equal(
    area_Polygons,
    area(SPolygonsDF)
  )    
  expect_equal(
    area_Lines,
    SpatialLinesLengths(SLinesDF)
  )  
})
