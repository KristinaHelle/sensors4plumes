#################################################################
# test plot.Simulations                                         #
#################################################################

test_that("plot.Simulations", {
  # data in memory
  data(radioactivePlumes_area)
  plot(radioactivePlumes_area, zcol = 3, main = "radioactive \n plumes \n area", col = bpy.colors())
  
  ## no numeric/factorial locations
  data(radioactivePlumes_local)
  plot(radioactivePlumes_local, zcol = 3)
  
  radioactivePlumes_noParam = radioactivePlumes_local
  radioactivePlumes_noParam@locations@data$index = NULL
  plot(radioactivePlumes_noParam, main = "no parameters")
  radioactivePlumes_manyParam = radioactivePlumes_area
  radioactivePlumes_manyParam@plumes = cbind(radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes,
                                             radioactivePlumes_area@plumes)
  radioactivePlumes_manyParam@locations@data = cbind(radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data,
                                                     radioactivePlumes_area@locations@data)
  plot(radioactivePlumes_manyParam, zcol = 3)
  
  # data not in memory
  bigSimulations = cbind(radioactivePlumes_area, 
                         radioactivePlumes_area, 
                         radioactivePlumes_area,
                         radioactivePlumes_area,
                         radioactivePlumes_area,
                         radioactivePlumes_area,
                         radioactivePlumes_area,
                         radioactivePlumes_area,
                         radioactivePlumes_area,
                         radioactivePlumes_area,
                         savePath = "bigSimulations",
                         overwrite = TRUE)
  #plot(bigSimulations, zcol = 3, maxpixels = 3e+7) # takes very long!
  
})

## for values integrate that clicking generated browsing through kinds (if desired, i.e. if several kinds given)