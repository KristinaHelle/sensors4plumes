#########################################################################
# test                        changeSimulationsPath                     #
#########################################################################

# inMemory -> nothing should happen
data(SimulationsSmall)
newSimulationsSmall = changeSimulationsPath(SimulationsSmall, "x.grd")

# not inMemory
## generate different raster files
x_12 = SimulationsSmall@values
x_21 = stack(SimulationsSmall@values[[2]], SimulationsSmall@values[[1]])
x_1 = SimulationsSmall@values[[1]]
x_2 = SimulationsSmall@values[[2]]
path1 = "/home/kristina/Desktop/s4p_article/test/"
writeRaster(x_1, file = paste0(path1, "/SimulationsSmall_a.grd"))
writeRaster(x_2, file = paste0(path1, "/SimulationsSmall_b.grd"))
writeRaster(x_12, file = paste0(path1, "/SimulationsSmall_ab.grd"))
writeRaster(x_21, file = paste0(path1, "/SimulationsSmall_ba.grd"))

plumes_i = data.frame(a = c("A", "A", "B", "B", "B"))
locations_i = SpatialIndexDataFrame(index = as.integer(1:9), data = data.frame(x = 1:9))
plumes_ii = rbind(plumes_i, plumes_i)
locations_ii = SpatialIndexDataFrame(index = as.integer(1:10), data = data.frame(x = 1:10))
values_ii_1 = raster(x = matrix(1:100, nrow = 10),
                crs = CRS("+init=epsg:4326"),
                xmn = -90, xmx = 90, ymn = -90, ymx = 90)
writeRaster(values_ii_1, file = paste0(path1, "/values_ii.grd"), overwrite = TRUE)
values_ii_2 = raster(paste0(path1, "/values_ii.grd")) # correct form for Simulations
values_ii_3 = raster(paste0(path.package("s4p"), # incorrect form for Simulations
                          "/extdata/fileFormats/raster/sourceA-date1_typea.tif"))
simulations_ii = Simulations(plumes = plumes_ii, locations = locations_ii, values = values_ii_2)
                             

# layer
simulations1a = Simulations(plumes = plumes_i, locations = locations_i, 
                            values =  raster(paste0(path1, "/SimulationsSmall_a.grd")))
simulations1b = Simulations(plumes = plumes_i, locations = locations_i, 
                            values =  raster(paste0(path1, "/SimulationsSmall_b.grd")))
simulations1B = changeSimulationsPath(simulations_i, path = paste0(path1, "/SimulationsSmall_b.grd"))

# brick
simulations1ab = Simulations(plumes = plumes_i, locations = locations_i, 
                             values =  brick(paste0(path1, "/SimulationsSmall_ab.grd"))) 
simulations1ba = Simulations(plumes = plumes_i, locations = locations_i, 
                             values =  brick(paste0(path1, "/SimulationsSmall_ba.grd"))) 
simulations1BA = changeSimulationsPath(simulations1ab, path = paste0(path1, "/SimulationsSmall_ba.grd"))

#stack
simulations2ab = Simulations(plumes = plumes_i, locations = locations_i, 
                             values = stack(c(paste0(path1, "/SimulationsSmall_a.grd"),
                                              paste0(path1, "/SimulationsSmall_b.grd")))) 
simulations2ba = Simulations(plumes = plumes_i, locations = locations_i, 
                             values = stack(c(paste0(path1, "/SimulationsSmall_b.grd"),
                                              paste0(path1, "/SimulationsSmall_a.grd")))) 
simulations2BA = changeSimulationsPath(simulations2ab, 
                                       path = c(paste0(path1, "/SimulationsSmall_b.grd"),
                                                paste0(path1, "/SimulationsSmall_a.grd")))

test_that("changeSimulationsPath",{
  
  expect_equal(SimulationsSmall, newSimulationsSmall) 
  
  expect_is(simulations1B@values, "RasterLayer")
  expect_equal(simulations1b, simulations1B)

  expect_is(simulations1BA@values, "RasterBrick")
  expect_equal(simulations1ba, simulations1BA)

  expect_is(simulations2BA@values, "RasterStack")
  expect_equal(simulations2ba, simulations2BA)

  # errors:
  ## file does not exist
  expect_error(changeSimulationsPath(simulations1a, path = ".grd"))
  ## wrong size
  data(RadioactivePlumes_area)
  expect_error(changeSimulationsPath(RadioactivePlumes_area, path = paste0(path1, "/SimulationsSmall_ab.grd")))   
  ## wrong format for use in Simulations

  expect_error(Simulations(locations = locations_ii,
                               plumes = plumes_ii,
                               values = values_ii))
  expect_error(changeSimulationsPath(simulations_ii,
                                     paste0(path.package("s4p"),
                                            "/extdata/fileFormats/raster/sourceA-date1_typea.tif")))
})


