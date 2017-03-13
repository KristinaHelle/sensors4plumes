##################################################################
# test        subsetSDF.SpatialPolygridDataFrame                 #
##################################################################

# combined parameters
SPolygridDFD1 = subsetSDF(SPolygridDFA1, grid_i = c(1,3,5), grid_j = c(2,4,6), 
                          coord_x = c(6,12), coord_y = c(2,12), locations = 5:12, data = "b")
# grid_ij
SPolygridDFD2 = subsetSDF(SPolygridDFD1, data = "b", grid_ij = matrix(c( TRUE, FALSE,  TRUE,
                                                                         FALSE, FALSE, FALSE, 
                                                                         FALSE, FALSE,  TRUE,
                                                                         FALSE, FALSE, FALSE,
                                                                         TRUE, FALSE,  TRUE),
                                                                      nrow = 5, byrow = TRUE))

# asymmetric input, boudaries through cells, margin cells empty (deleted)
SPolygridDFD3 = subsetSDF(SPolygridDFB2, coord_x = c(7,12), coord_y = c(5,11)) 

# point result 
SPolygridDFD4 = subsetSDF(SPolygridDFD2, coord_x = c(9,12), coord_y = c(3,5)) 

test_that("subsetSDF.SpatialPolygridDataFrame", {
  expect_equal(
    SPolygridDFD2, SPolygridDFD1
  )
  
  expect_equal(
    dename(bbox(SPolygridDFD2), kind = "dimnames"), 
    matrix(c(6,2,12,12), nrow = 2)
  )
  
#  expect_is(
#    SPolygridDFD4, 
#    "SpatialPointsDataFrame"
#  )
  
  # parameters to be ignored (-> warnings)
  expect_warning(
    subsetSDF(SPolygridDFA1, data = "b", 
              grid_i = c(1,3,5,6), grid_j = c(1,3,4), 
              coord_x = c(-2,11), coord_y = c(1,13), index = 1:7,
              grid_ij = matrix(c( TRUE, FALSE,  TRUE,
                                  FALSE, FALSE, FALSE, 
                                  FALSE, FALSE,  TRUE,
                                  FALSE, FALSE, FALSE,
                                  TRUE, FALSE,  TRUE),
                               nrow = 5, byrow = TRUE))
  )
  # result linear (-> SpatialPointsDataFrame, warning)
#  expect_warning(
#    subsetSDF(SPolygridDFA1, coord_x = c(10,12), coord_y = c(3,11)) 
#  )
  # result empty (-> error)
  expect_error(
    subsetSDF(SPolygridDFA1, coord_x = c(8,12), coord_y = c(3,11), locations = c(5,6,12)) 
  )
})



