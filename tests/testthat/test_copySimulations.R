#########################################################################
# test                        copySimulations                           #
#########################################################################

simulations_ii # from test_changeSimulationsPath
dir.create(paste0(path1, "/test1"))

save(simulations_ii, file = paste0(path1, "/test1/simulations_ii.Rdata"))

RadioactivePlumes_area2 = 
  copySimulations(RadioactivePlumes_area, newPath = paste0(path1, "test1"))

simulations_iii = copySimulations(simulations_ii, newFile = "simulations_ii")

test_that("copySimulations", {
  # values are kept
  i = sample.int(nlayers(RadioactivePlumes_area2@values),1)
  j = sample.int(nrow(RadioactivePlumes_area2@values[[1]]),1)
  expect_equal(
    getValues(RadioactivePlumes_area2@values[[i]],row = j, nrows = 1),  
    getValues(RadioactivePlumes_area@values[[i]],row = j, nrows = 1)  
  )
  # paths are changed
  fileName = rev(unlist(strsplit(RadioactivePlumes_area2@values[[i]]@file@name, "/")))[[1]]
  expect_equal( # same file name
    fileName,
    rev(unlist(strsplit(RadioactivePlumes_area@values[[i]]@file@name, "/")))[[1]]
  )
  expect_equal( # new path
    strsplit(RadioactivePlumes_area2@values[[i]]@file@name, split = paste0("/", fileName))[[1]],
    paste0(path1, "test1")
  )
  
  expect_true(
    file.exists("simulations_ii.Rdata")
  )  
  
  expect_error(
    copySimulations(simulations_ii, newPath = "simulations_ii")  
  )
})
