################################################################
# test              points2polygrid                            #
################################################################
#source('tests/testthat/test_SpatialPolygridDataFrame_prepare.R')
test_that("points2polygrid point", {
  Grid = coordinatesA1
  roundReg = runif(24, -0.0025, 0.0025)
  irreg = runif(24)
  testPoints = list()
  testPoints[["Grid00"]] = Grid
  testPoints[["Grid11"]] = Grid + roundReg
  testPoints[["Grid22"]] = Grid + irreg
  testPoints[["Grid10"]] = cbind(Grid[,1] + roundReg[1:12], Grid[,2])
  testPoints[["Grid01"]] = cbind(Grid[,1],                  Grid[,2] + roundReg[13:24])
  testPoints[["Grid20"]] = cbind(Grid[,1] + irreg[1:12],    Grid[,2])
  testPoints[["Grid02"]] = cbind(Grid[,1],                  Grid[,2] + irreg[13:24])
  testPoints[["Grid12"]] = cbind(Grid[,1] + roundReg[1:12], Grid[,2] + irreg[13:24])
  testPoints[["Grid21"]] = cbind(Grid[,1] + irreg[1:12],    Grid[,2] + roundReg[13:24])
  testPoints[["LineX00"]] = cbind(Grid[,1], 0)
  testPoints[["LineX11"]] = cbind(Grid[,1], 0) + roundReg
  testPoints[["LineX22"]] = cbind(Grid[,1], 0) + irreg
  testPoints[["LineX10"]] = cbind(Grid[,1] + roundReg[1:12], 0)
  testPoints[["LineX01"]] = cbind(Grid[,1],                  0 + roundReg[13:24])
  testPoints[["LineX20"]] = cbind(Grid[,1] + irreg[1:12],    0)
  testPoints[["LineX02"]] = cbind(Grid[,1],                  0 + irreg[13:24])
  testPoints[["LineX12"]] = cbind(Grid[,1] + roundReg[1:12], 0 + irreg[13:24])
  testPoints[["LineX21"]] = cbind(Grid[,1] + irreg[1:12],    0 + roundReg[13:24])
  testPoints[["LineY00"]] = testPoints[["LineX00"]][,2:1]
  testPoints[["LineY11"]] = testPoints[["LineX11"]][,2:1]
  testPoints[["LineY22"]] = testPoints[["LineX22"]][,2:1]
  testPoints[["LineY10"]] = testPoints[["LineX10"]][,2:1]
  testPoints[["LineY01"]] = testPoints[["LineX01"]][,2:1]
  testPoints[["LineY20"]] = testPoints[["LineX20"]][,2:1]
  testPoints[["LineY02"]] = testPoints[["LineX02"]][,2:1]
  testPoints[["LineY12"]] = testPoints[["LineX12"]][,2:1]
  testPoints[["LineY21"]] = testPoints[["LineX21"]][,2:1]
  testPoints[["Point00"]] = cbind(rep(0,12), 0)
  testPoints[["Point11"]] = cbind(rep(0,12), 0) + roundReg
  testPoints[["Point22"]] = cbind(rep(0,12), 0) + irreg
  testPoints[["Point10"]] = cbind(0 + roundReg[1:12], 0)
  testPoints[["Point01"]] = cbind(0,                  0 + roundReg[13:24])
  testPoints[["Point20"]] = cbind(0 + irreg[1:12],    0)
  testPoints[["Point02"]] = cbind(0,                  0 + irreg[13:24])
  testPoints[["Point12"]] = cbind(0 + roundReg[1:12], 0 + irreg[13:24])
  testPoints[["Point21"]] = cbind(0 + irreg[1:12],    0 + roundReg[13:24])          
  # 2dim, 1dim, 0dim grid before/after rounding
  testResult = list()
  for (i in c(1:9, 12, 16, 17, 21, 25, 26, 30)){
    testResult[[names(testPoints)[i]]] = points2polygrid(points = testPoints[[i]], tolerance = 0.005)
  }
  for (i in c(10, 11, 13:15, 18:20, 22:24, 27, 33:36)){
    expect_warning(
      testResult[[names(testPoints)[i]]] <- points2polygrid(points = testPoints[[i]], tolerance = 0.005)
    )
  }
  for (i in 28:29){
    expect_error(
      testResult[[names(testPoints)[i]]] <- try(points2polygrid(points = testPoints[[i]], tolerance = 0.005))
    ) 
  }  
  for (i in 31:32){
    expect_warning(
      expect_error(
        testResult[[names(testPoints)[i]]] <- try(points2polygrid(points = testPoints[[i]], tolerance = 0.005))
      ) 
    )
  }
})

# points given: regular, regular telescopic, irregular; wrong (wrong kind, too few columns; on a line)
# points & grid given: irregular, non-rectangular 
#  (proj: not, only gr, only p, both contradictory); (overlapping: total, partially, not at all)
# grid & index given (length of index fits / does not fit)
# no data, correct data, data of wrong length
# varying tolerance (irregular points, no grid)

test_that("points2polygrid", {  
  ## wrong points
  expect_warning(
    points2polygrid(# under this tolerance, points are in one line
      points = coordinates2[1:6,],
      tolerance = 2)
  )
  
  expect_warning(
    points2polygrid(# multiple points with exactly identical coordinates
      points = rbind(coordinates2, coordinates2))
    )
  )
  
  expect_error(
    points2polygrid(# only one column
      points = matrix(coordinates1, ncol = 1))
  )
  
  expect_error(
    points2polygrid(# non-numeric 
      points = data.frame(x1 = c("x", "y", "z"), x2 = c("A", "B", "C"))))
  
  
  ## regular
  SPolygridDFF1 = points2polygrid(
    points = coordinates1,
    data = dataFrame3)
  
  expect_equal(# regular points become grid coordinates
    coordinates1,
    coordinates(SPolygridDFF1))
  
  SPolygridDFF2 = points2polygrid(
    points = coordinates1[36:1,],
    data = dataFrame3)
  
  expect_equal(# order of points matters
    SPolygridDFF1@data[SPolygridDFF1@index,],
    SPolygridDFF2@data[SPolygridDFF2@index,][36:1,])
  
  
  ## rectangular
  SPolygridDFF3 = points2polygrid(
    points = coordinates3)
  
  expect_equal(
    SPolygridDFF3@grid,
    grid4)
  
  ## telescopic
  SPolygridDFF4 = points2polygrid(
    points = coordinatesA1)
  
  expect_equal(# regular -> resolution equals smallest distance between points
    SPolygridDFF4@grid@cellsize,
    c(1,1))
  
  SPolygridDFF5 = points2polygrid(
    points = coordinates4)
  
  ## irregular
  SPolygridDFF6 = points2polygrid(
    points = coordinates2,
    data = dataFrame3
  )
  
  expect_less_than(
    sum(SPolygridDFF6@grid@cells.dim),
    2500)
  
  
  # points and grid
  expect_warning(
    SPolygridDFF7 <- points2polygrid(# conflicting projections 
      points = spatialPoints1,
      grid = spatialGrid2))
  
  expect_equal(
    SPolygridDFF7@proj4string,
    spatialGrid1@proj4string)
  
  SPolygridDFF8 = points2polygrid(# irregular points in grid
    points = coordinates2,
    grid = spatialGrid1)
  
  expect_equal(
    SPolygridDFF8@grid,
    spatialGrid1@grid)
  
  expect_warning(
    SPolygridDFF9 <- points2polygrid(# grid non-overlapping with points, merge all, data of closest point
      points = coordinates1,
      grid = grid3,
      data = dataFrame3))
  
  expect_equal(
    SPolygridDFF9@data,
    dataFrame3[6, ])
  
  # points, grid & index -> points ignored
  expect_warning(
    SPolygridDFF10 <- points2polygrid(
      points = coordinatesA1[12:1,],
      grid = grid1,
      index = index1,
      data = dataFrame1)
  )
  
  expect_equal(
    SPolygridDFF10,
    SPolygridDFA1
  )
  
  # points & data, wrong number
  expect_warning(
    SPolygridDFF11 <- points2polygrid(
      points = coordinatesA1,
      data = dataFrame3)
  )
  expect_equal(
    SPolygridDFF11@data,
    dataFrame3[1:12,]
  )
  
  # index & data, index not 1:nrow(data)
  expect_warning(
    SPolygridDFF12 <- points2polygrid(
      grid = grid1,
      index = as.integer(index1 + 6),
      data = dataFrame1)  
  )
  
  expect_equal(
    SPolygridDFF12@data,
    dataFrame1[7:12,]
  )
  
  expect_warning(
    SPolygridDFF13 <- points2polygrid(
      grid = grid1,
      index = as.integer(index1 - 6),
      data = dataFrame1)  
  )
  expect_equal(
    SPolygridDFF13@data,
    dataFrame1[1:6,]
  )
  
  indexF1 = index1
  indexF1[index1 <= 6] = NA
  indexF2 = index1
  indexF2[index1 == 4] = as.integer(-4)
  expect_warning(
    SPolygridDFF14 <- points2polygrid(
      grid = grid1,
      index = indexF1,
      data = dataFrame1)  
  )
  expect_equal(
    SPolygridDFF14@index,
    SPolygridDFF13@index
  )
  expect_warning(
    SPolygridDFF15 <- points2polygrid(
      grid = grid1,
      index = indexF2,
      data = dataFrame1)  
  )
  expect_equal(
    SPolygridDFF15@data,
    dataFrame1[-4,])
  
  expect_warning(
    SPolygridDFF16 <- points2polygrid(
      grid = grid1,
      index = as.integer(1:36),
      data = dataFrame1)  
  )
  expect_equal(
    SPolygridDFF16@index[13:36],
    as.numeric(rep(NA, 24))
  )
  expect_equal(
    SPolygridDFF16@data,
    dataFrame1
  )
  expect_warning(
    expect_error(
      points2polygrid(
        grid = grid1,
        index = as.integer(1:12),
        data = dataFrame1)  
    )
  )
  expect_warning(
    expect_error(
      points2polygrid(
        grid = grid1,
        index = as.integer(rep(1:12, 4)),
        data = dataFrame1)  
    )  
  )



})
