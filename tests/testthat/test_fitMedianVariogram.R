###############################################################
# test fitMedianVariogram                                     #
###############################################################
# if simulations@locations is of different class
# is locations are subsetted or not

test_that("fitMedianVariogram", {
  data(radioactivePlumes_local)
  # samples
  samplePlumes = sample.int(nPlumes(radioactivePlumes_local), 10)
  sampleLocations = sample.int(nLocations(radioactivePlumes_local), 500)
  
  # tun locations into SDF of different class
  # SpatialIndexDataFrame
  SIndexDF1 = SpatialIndexDataFrame(index = 1:nLocations(radioactivePlumes_local),
                                    data = radioactivePlumes_local@locations@data)
  radioactivePlumes_index = radioactivePlumes_local
  radioactivePlumes_index@locations = SIndexDF1
  expect_error(
    vgm_index <- fitMedianVariogram(simulations = radioactivePlumes_index,
                                   plumes  = samplePlumes,
                                   values = 1)  
  )
  # SpatialPolygridDataFrame
  radioactivePlumes_polygrid  = subset(radioactivePlumes_local, locations = 1:2001)
  vgm_polygridAll <- fitMedianVariogram(simulations = radioactivePlumes_polygrid,
                                  plumes  = samplePlumes,
                                  values = 1) 
  vgm_polygridSubset <- fitMedianVariogram(simulations = radioactivePlumes_polygrid,
                                           plumes  = samplePlumes,
                                           locations = sampleLocations,
                                           values = 1) 
})