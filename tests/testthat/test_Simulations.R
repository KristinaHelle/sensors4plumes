###################################################################
# test Simulations()                                              #
###################################################################
## Simulations()
# wrong format of input
# size does not fit: locations, plumes
# small (all kind of SDF); big (raster from file: brick, stack)

## nLocations, nPlumes

# load data: can be removed if test_dataArtificial is run before
data(SIndexDF)
data(SPointsDF)
data(SPixelsDF)
data(SPolygridDF)
data(SPolygonsDF)
data(SLinesDF)

# generate SGridDF (not SDF!)
SGridDF = SPixelsDF
fullgrid(SGridDF) = TRUE

# generate plumes and values
plumes1 = data.frame(source = c("A", "A", "B", "B", "B"),
                     date = c("2000-01-01", "2000-04-01", "2000-07-01", "2000-01-01", "2000-01-03"),
                     totalCost = runif(5, min = 5, max = 15))

values1 = stack(raster(x = replicate(n = 5, expr = rlnorm(2, sdlog = 2)), 
                      xmn = -90, xmx = 90, ymn = -90, ymx = 90, 
                      crs = "+init=epsg:4326" ), 
               raster(x = replicate(n = 5, expr = rnorm(2, m = 10, sd = 3)), 
                      xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                      crs = "+init=epsg:4326" ))
values2 = stack(raster(x = t(matrix(0:9, nrow = 5)),
                xmn = -90, xmx = 90, ymn = -90, ymx = 90, 
                crs = "+init=epsg:4326"),
                raster(x = t(matrix(seq(10, 100, 10), nrow = 5)), 
                       xmn = -90, xmx = 90, ymn = -90, ymx = 90, 
                       crs = "+init=epsg:4326"))

# simulations
SimulationsIndex1 = Simulations(locations = subsetSDF(SIndexDF, locations = 1:2),
                               plumes = plumes1,
                               values = values1)
SimulationsIndex2 = Simulations(locations = subsetSDF(SIndexDF, locations = 1:2),
                                plumes = plumes1,
                                values = values2)
SimulationsPoints1 = Simulations(locations = subsetSDF(SPointsDF, locations = 1:2),
                                plumes = plumes1,
                                values = values1)
SimulationsPixels1 = Simulations(locations = subsetSDF(SPixelsDF, locations = 1:2),
                                plumes = plumes1,
                                values = values1)
SimulationsPolygrid1 = Simulations(locations = subsetSDF(SPolygridDF, locations = 1:2),
                                plumes = plumes1,
                                values = values1)
SimulationsPolygons1 = Simulations(locations = subsetSDF(SPolygonsDF, locations = 1:2),
                                plumes = plumes1,
                                values = values1)
SimulationsLines1 = Simulations(locations = subsetSDF(SLinesDF, locations = 1:2),
                                plumes = plumes1,
                                values = values1)
test_that("Simulations()", {
  expect_is(
    SimulationsIndex1, "Simulations"
  )
  expect_is(
    SimulationsPoints1, "Simulations"
  )
  expect_is(
    SimulationsPixels1, "Simulations"
  )
  expect_is(
    SimulationsPolygrid1, "Simulations"
  )
  expect_is(
    SimulationsPolygons1, "Simulations"
  )
  expect_is(
    SimulationsLines1, "Simulations"
  )
  
  expect_error( # wrong class of locations (not SDF): 
    Simulations(locations = SGridDF[1:2,],
                plumes = plumes1, values = values1)
  )  
  
  expect_error( # wrong size: locations vs. values
    Simulations(locations = SIndexDF,
                plumes = plumes1, values = values1)
  )
  
  expect_error( # wrong size: plumes vs. values
    Simulations(locations = subsetSDF(SIndexDF, locations = 1:2),
                plumes = plumes1[1:4,], values = values1)
  )
})


# ----------------- nLocations & nPlumes & nKinds ------------------------ #
data(SimulationsSmall) # replace by final data

nLocSS = nLocations(SimulationsSmall)
nPlSS = nPlumes(SimulationsSmall)
nKSS = nKinds(SimulationsSmall)
test_that("nLocations.Simulations, nPlumes.Simulations, nKinds.Simulations", {
  expect_equal(
    nLocSS, nrow(SimulationsSmall@values)
  )  
  expect_equal(
    nPlSS, ncol(SimulationsSmall@values)
  ) 
  expect_equal(
    nKSS, nlayers(SimulationsSmall@values)
  )
})

