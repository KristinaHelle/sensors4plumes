library(s4p)
library(testthat)


fun1 = idw0z #replaceDefault(idw0, newDefaults = list(formula = z ~ 1), type = "fun_interpolation.interpolate")[[1]] 
fun2 = replaceDefault(krige0, newDefaults = list(formula = z ~ 1, beta = NA, model = vgm(10, "Exp", 50000, 5), ... = NA), type = "fun_interpolation.interpolate")[[1]]
fun3 = replaceDefault(krige0, newDefaults = list(formula = z ~ 1, beta = NA, model = vgm(100, "Sph", 500, 50), ... = NA), type = "fun_interpolation.interpolate")[[1]]


# test with different SDF
#data(SIndexDF) # does not work, SIndexDF does not have spatial properties
data(SPointsDF)
data(SPixelsDF)
#data(SPolygridDF) # does not work directly as krige0 / idw0 do only work on spatial objects; can be solved by transforming polygrid into SPoints beforehand
data(SLinesDF)
data(SPolygonsDF)

# subsetting
dS1_0 = 7:9
dS2_0 = 1:6
nS1_0 = 1:3
nS2_0 = 4:9
d_0 = c(1,2,3,8,9)
n_0 = 2:9

# Points
intPoints_1 = interpol0reg(fun_interpolation = fun2,
                            dataSplit = list(dS1_0, dS2_0),
                            newdataSplit = list(nS1_0, nS2_0),
                            y = as.matrix(SPointsDF@data)[d_0,, drop = FALSE], 
                            data = SPointsDF[d_0,], 
                            newdata = SPointsDF[n_0,],
                            dataLoc = d_0, 
                            newdataLoc = n_0)

dimnames(intPoints_1)[[1]] = n_0

# -  test - 
intPoints_part1 = fun2(y = as.matrix(SPointsDF@data)[intersect(dS1_0, d_0),, drop = FALSE], 
                       data = SPointsDF[intersect(dS1_0, d_0),],
                       newdata = SPointsDF[n_0,])
intPoints_part2 = fun2(y = as.matrix(SPointsDF@data)[intersect(dS2_0, d_0),, drop = FALSE], 
                       data = SPointsDF[intersect(dS2_0, d_0),],
                       newdata = SPointsDF[n_0,])
intPoints_parts = intPoints_part1
intPoints_parts[is.element(n_0, nS2_0),] = intPoints_part2[is.element(n_0, nS2_0),]

expect_equivalent(intPoints_1, intPoints_parts)

# Pixels
intPixels_1 = interpol0reg(fun_interpolation = list(fun1, fun2),
                                  dataSplit = list(dS1_0, dS2_0),
                                  newdataSplit = list(nS1_0, nS2_0),
                                  y = as.matrix(SPixelsDF@data)[d_0,, drop = FALSE], 
                                  data = SPixelsDF[d_0,], 
                                  newdata = SPixelsDF[n_0,], 
                                  dataLoc = d_0, 
                                  newdataLoc = n_0)

intPixels_part1 = fun1(y = as.matrix(SPixelsDF@data)[intersect(dS1_0, d_0),, drop = FALSE], 
                       data = SPixelsDF[intersect(dS1_0, d_0),],
                       newdata = SPixelsDF[intersect(nS1_0, n_0),])
intPixels_part2 = fun2(y = as.matrix(SPixelsDF@data)[intersect(dS2_0, d_0),, drop = FALSE], 
                       data = SPixelsDF[intersect(dS2_0, d_0),],
                       newdata = SPixelsDF[intersect(nS2_0, n_0),])
intPixels_parts = matrix(nrow = nrow(intPixels_1), ncol(intPixels_1))
intPixels_parts[is.element(n_0, nS1_0),] = intPixels_part1
intPixels_parts[is.element(n_0, nS2_0),] = intPixels_part2

expect_equivalent(intPixels_1, intPixels_parts)

# Lines: idw0 fails
#idw0(formula = z ~ 1, y = as.matrix(SLinesDF2@data)[1:2,], data = SLinesDF2[1:2,],  newdata = SLinesDF2) # Error (=> not MY problem)
# Polygons: idw0 fails
#idw0(formula = z ~ 1, y = as.matrix(SPolygonsDF@data)[1:2,, drop = FALSE], data = SPolygonsDF[1:2],  newdata = SPolygonsDF) # Error (=> not MY problem)


# newdata not completely covered and overlapping
dS1_nc0 = 7:9
dS2_nc0 = 1:6
nS1_nc0 = 1:3
nS2_nc0 = c(2,5:9)
d_nc0 = c(1,2,3,8,9)
n_nc0 = 1:8
expect_warning(# waring about overlapping newdataSplit
  int_nc1 <- interpol0reg(fun_interpolation = fun2,
                         dataSplit = list(dS1_nc0, dS2_nc0),
                         newdataSplit = list(nS1_nc0, nS2_nc0),
                         y = as.matrix(SPointsDF@data)[d_nc0,, drop = FALSE], 
                         data = SPointsDF[d_nc0,], 
                         newdata = SPointsDF[n_nc0,], 
                         dataLoc = d_nc0, 
                         newdataLoc = n_nc0)  
)
dimnames(int_nc1)[[1]] = n_nc0

expect_equivalent(int_nc1[4], as.integer(NA)) # not in any nSX
expect_false(int_nc1[2] == intPoints_1[1]) # in both nSX, determined by last interpolation
expect_equivalent(int_nc1[c(3, 5:8),], intPoints_1[c(2, 4:7),]) # other values are the same 

# reversed, multiple, invalid, NA indices
# reversed
dS1_r01 = c(9,7,8)
dS2_r01 = c(3,2,6,5,1,4)
nS1_r01 = c(3,1,2)
nS2_r01 = c(7,5,8,9,6,4)
d_r01 = c(1,9,3,8,2)
n_r01 = c(3,4,7,6,8,9,2,5)

int_r01 = interpol0reg(fun_interpolation = fun2,
                                  dataSplit = list(dS1_r01, dS2_r01),
                                  newdataSplit = list(nS1_r01, nS2_r01),
                                  y = as.matrix(SPointsDF@data)[d_r01,, drop = FALSE], 
                                  data = SPointsDF[d_r01,], 
                                  newdata = SPointsDF[n_r01,], 
                                  dataLoc = d_r01, 
                                  newdataLoc = n_r01)
dimnames(int_r01)[[1]] = n_r01

expect_equivalent(int_r01[order(dimnames(int_r01)[[1]]),,drop = FALSE], intPoints_1) # changing order of input does not matter; output order as required

# multiple
dS1_m01 = c(7:9,9)
dS2_m01 = c(1,1:6)
nS1_m01 = c(1:3,3)
nS2_m01 = c(4,4:8)
d_m01 = c(1,1,2,3,8,8,9)#c(1,1,2,3,4,8,8,9)
n_m01 = c(2,2:9,9)

expect_error(# this example does not work (multiple input is not the problem)
  fun2(y = as.matrix(SPointsDF@data)[c(1,2,2,3,4,8,8,9)[is.element(c(1,2,2,3,4,8,8,9), dS2_m01)],, drop = FALSE], 
       data = SPointsDF[c(1,2,2,3,4,8,8,9)[is.element(c(1,2,2,3,4,8,8,9), dS2_m01)],], 
       newdata = SPointsDF[n_m01[is.element(n_m01, nS2_m01)],])
)

int_m01 = interpol0reg(fun_interpolation = fun2,
                              dataSplit = list(dS1_m01, dS2_m01),
                              newdataSplit = list(nS1_m01, nS2_m01),
                              y = as.matrix(SPointsDF@data)[d_m01,, drop = FALSE], 
                              data = SPointsDF[d_m01,], 
                              newdata = SPointsDF[n_m01,], 
                              dataLoc = d_m01, 
                              newdataLoc = n_m01)
dimnames(int_m01)[[1]] = n_m01

expect_equivalent(int_m01[2:8], intPoints_1[1:7]) # multiple input has no effect
expect_equivalent(int_m01[1,],int_m01[2,]) # multiple output has same value 
expect_equivalent(int_m01[9,],int_m01[10,]) 

# containing NA: throws errors when used to subset Spatial objects
dS1_n01 = c(NA, 7:9)
dS2_n01 = c(1:6, NA)
nS1_n01 = c(1:3, NA, NA)
nS2_n01 = c(NA, 4:9)
d_n01 = c(NA, NA, 1,2,3,8,9)
n_n01 = c(2:9, NA, NA)

expect_error(# input with NA: NA not permitted in row index
  fun2(y = as.matrix(SPointsDF@data)[d_n01[is.element(d_n01, dS1_n01)],, drop = FALSE], 
       data = SPointsDF[d_n01[is.element(d_n01, dS1_n01)],], 
       newdata = SPointsDF[n_0[is.element(n_0, nS1_0)],])
)
expect_error(# output with NA: NA not permitted in row index
  fun2(y = as.matrix(SPointsDF@data)[d_0[is.element(d_0, dS1_0)],, drop = FALSE], 
       data = SPointsDF[d_0[is.element(d_0, dS1_0)],], 
       newdata = SPointsDF[n_n01[is.element(n_n01, nS2_n01)],])
)

# containing invalid indices
dS1_i01 = c(0, 7:9, 10)
dS2_i01 = c(-1, 1:6)
nS1_i01 = c(11, 1:3)
nS2_i01 = c(-3, 4:9)
d_i01 = c(-1, 0, 1,2,3,8,9, 10)
n_i01 = c(2:9, -3, 11)

expect_error(# input with invalid indices: subscript out of bounds
  fun2(y = as.matrix(SPointsDF@data)[d_i01[is.element(d_i01, dS1_i01)],, drop = FALSE], 
       data = SPointsDF[d_i01[is.element(d_i01, dS1_i01)],], 
       newdata = SPointsDF[n_0[is.element(n_0, nS1_0)],])
)
expect_error(# output with invalid indices: only 0's may be mixed with negative subscripts
  fun2(y = as.matrix(SPointsDF@data)[d_0[is.element(d_0, dS1_0)],, drop = FALSE], 
       data = SPointsDF[d_0[is.element(d_0, dS1_0)],], 
       newdata = SPointsDF[n_i01[is.element(n_i01, nS2_i01)],])
)

# some parameters given as lists (of equal elements), others not
expect_warning(
  intPoints_lll <- interpol0reg(fun_interpolation = list(fun2, fun2),
                                      dataSplit = list(dS1_0, dS1_0),
                                      newdataSplit = list(nS1_0, nS1_0),
                                      y = as.matrix(SPointsDF@data)[d_0,, drop = FALSE], 
                                      data = SPointsDF[d_0,], 
                                      newdata = SPointsDF[n_0,], 
                                      dataLoc = d_0, 
                                      newdataLoc = n_0)  
)
expect_warning(
  intPoints_ell <- interpol0reg(fun_interpolation = fun2,
                                       dataSplit = list(dS1_0, dS1_0),
                                       newdataSplit = list(nS1_0, nS1_0),
                                       y = as.matrix(SPointsDF@data)[d_0,, drop = FALSE], 
                                       data = SPointsDF[d_0,], 
                                       newdata = SPointsDF[n_0,], 
                                       dataLoc = d_0, 
                                       newdataLoc = n_0)  
)
expect_warning(
  intPoints_lel <- interpol0reg(fun_interpolation = list(fun2, fun2),
                                       dataSplit = dS1_0,
                                       newdataSplit = list(nS1_0, nS1_0),
                                       y = as.matrix(SPointsDF@data)[d_0,, drop = FALSE], 
                                       data = SPointsDF[d_0,], 
                                       newdata = SPointsDF[n_0,], 
                                       dataLoc = d_0, 
                                       newdataLoc = n_0)  
)
intPoints_lle <- interpol0reg(fun_interpolation = list(fun2, fun2),
                                       dataSplit = list(dS1_0, dS1_0),
                                       newdataSplit = nS1_0,
                                       y = as.matrix(SPointsDF@data)[d_0,, drop = FALSE], 
                                       data = SPointsDF[d_0,], 
                                       newdata = SPointsDF[n_0,], 
                                       dataLoc = d_0, 
                                       newdataLoc = n_0) 
intPoints_eee <- interpol0reg(fun_interpolation = fun2,
                                     dataSplit = dS1_0,
                                     newdataSplit = nS1_0,
                                     y = as.matrix(SPointsDF@data)[d_0,, drop = FALSE], 
                                     data = SPointsDF[d_0,], 
                                     newdata = SPointsDF[n_0,], 
                                     dataLoc = d_0, 
                                     newdataLoc = n_0)  
expect_equivalent(intPoints_lll, intPoints_ell)
expect_equivalent(intPoints_lll, intPoints_lel)
expect_equivalent(intPoints_lll, intPoints_lle)
expect_equivalent(intPoints_lll, intPoints_eee)

expect_warning(
  expect_error(# lists of different lengths
    interpol0reg(fun_interpolation = list(fun2),
                        dataSplit = list(dS1_0, dS1_0),
                        newdataSplit = list(nS1_0, nS1_0),
                        y = as.matrix(SPointsDF@data)[d_0,, drop = FALSE], 
                        data = SPointsDF[d_0,], 
                        newdata = SPointsDF[n_0,], 
                        dataLoc = d_0, 
                        newdataLoc = n_0)  
  )
)
expect_warning(
  expect_error(# lists of different lengths
    interpol0reg(fun_interpolation = list(fun2),
                        dataSplit = list(dS1_0, dS1_0),
                        newdataSplit = list(nS1_0, nS1_0),
                        y = as.matrix(SPointsDF@data)[d_0,, drop = FALSE], 
                        data = SPointsDF[d_0,], 
                        newdata = SPointsDF[n_0,], 
                        dataLoc = d_0, 
                        newdataLoc = n_0)  
  )
)
  expect_error(# lists of different lengths
    interpol0reg(fun_interpolation = list(fun2, fun2),
                        dataSplit = list(dS1_0, dS2_0),
                        newdataSplit = list(nS1_0),
                        y = as.matrix(SPointsDF@data)[d_0,, drop = FALSE], 
                        data = SPointsDF[d_0,], 
                        newdata = SPointsDF[n_0,], 
                        dataLoc = d_0, 
                        newdataLoc = n_0)  
  )

# test for ERRORS
# parameters of wrong classes
expect_error( # fun is no function
  interpol0reg(fun_interpolation = dS2_0,
                      dataSplit = fun2,
                      newdataSplit = nS1_0,
                      y = as.matrix(SPointsDF@data)[d_0,, drop = FALSE], 
                      data = SPointsDF[d_0,], 
                      newdata = SPointsDF[n_0,], 
                      dataLoc = d_0, 
                      newdataLoc = n_0) 
)



# lengths of y, data, dataLoc are not the same
expect_error( # y has wrong length
  interpol0reg(fun_interpolation = fun2,
                      dataSplit = dS2_0,
                      newdataSplit = nS1_0,
                      y = as.matrix(SPointsDF@data)[d_0[1:4],, drop = FALSE], 
                      data = SPointsDF[d_0,], 
                      newdata = SPointsDF[n_0,], 
                      dataLoc = d_0, 
                      newdataLoc = n_0) 
)
expect_error( # data has wrong length
  interpol0reg(fun_interpolation = fun2,
                      dataSplit = dS2_0,
                      newdataSplit = nS1_0,
                      y = as.matrix(SPointsDF@data)[d_0,, drop = FALSE], 
                      data = SPointsDF[d_0[1:3],], 
                      newdata = SPointsDF[n_0,], 
                      dataLoc = d_0, 
                      newdataLoc = n_0) 
)
expect_error( # dataLoc has wrong length
  interpol0reg(fun_interpolation = fun2,
                      dataSplit = dS2_0,
                      newdataSplit = nS1_0,
                      y = as.matrix(SPointsDF@data)[d_0,, drop = FALSE], 
                      data = SPointsDF[d_0,], 
                      newdata = SPointsDF[n_0,], 
                      dataLoc = d_0[2:5], 
                      newdataLoc = n_0) 
)
# lengths of newdata, newdataLoc are not the same
expect_error( # newdata, newdataLoc have different length
  interpol0reg(fun_interpolation = fun2,
                      dataSplit = dS2_0,
                      newdataSplit = nS1_0,
                      y = as.matrix(SPointsDF@data)[d_0,, drop = FALSE], 
                      data = SPointsDF[d_0,], 
                      newdata = SPointsDF[n_0,], 
                      dataLoc = d_0, 
                      newdataLoc = n_0[1:7]) 
)




# test with bigger (realistic) example
data(radioactivePlumes_local)
radPlumes_small = subsetSDF(extractSpatialDataFrame(radioactivePlumes_local, plumes = 1:20, kinds = 1), locations = 1:2001)
# locations must be transformed from Polygrid to Pixels or Points
radPlumes_smallP = as(radPlumes_small, "SpatialPointsDataFrame")
names(radPlumes_smallP)[1] = "z"

locIn1 = as.integer(seq(1,2001, 30))
y1 = as.matrix(radPlumes_small@data[locIn1,])
data1 = radPlumes_smallP[locIn1,]
newdata1 = radPlumes_smallP

fun1 = idw0z #replaceDefault(idw0, newDefaults = list(formula = z ~ 1), type = "fun_interpolation.interpolate")[[1]] 
fun2 = replaceDefault(krige0, newDefaults = list(formula = z ~ 1, beta = NA, model = vgm(10, "Exp", 50000, 5), ... = NA), type = "fun_interpolation.interpolate")[[1]]
fun3 = replaceDefault(krige0, newDefaults = list(formula = z ~ 1, beta = NA, model = vgm(100, "Sph", 500, 50), ... = NA), type = "fun_interpolation.interpolate")[[1]]


# test: missing parameters
# parameters with default values: fun, dataSplit, newdataSplit, newdataLoc
int__fun = interpol0reg(
                      y = y1, 
                      data = data1, 
                      newdata = newdata1, 
                      dataLoc = locIn1)  
int__fun1 = interpol0reg(fun = fun1,
                                dataSplit = 1:2001,
                                newdataSplit = 1:2001,
                                y = y1, 
                                data = data1, 
                                newdata = newdata1, 
                                dataLoc = locIn1, 
                                newdataLoc = 1:2001) 
expect_equivalent(int__fun, int__fun1)


# test with true splitting
dS1_1 = 1:500
dS2_1 = 501:1000
dS3_1 = 1001:2001
nS1_1 = 1:1001
nS2_1 = 1:1501
nS3_1 = 501:2001
d_1 = seq(1, 2001, 25)
n_1 = 1:1501  

interpol1 = replaceDefault(interpol0reg,
                           newDefaults = list(
                             fun_interpolation = list(fun1, fun2, fun3),
                             dataSplit = list(dS1_1, dS2_1, dS3_1),
                             newdataSplit = list(nS1_1, nS2_1, nS3_1)), 
                           type = "fun_interpolationSplit.interpolate")

expect_warning( # overlapping newdataSplit
  int1 <- interpol1[[1]](y = as.matrix(radPlumes_small@data[d_1,]), 
                        data = radPlumes_smallP[d_1,], 
                        newdata = radPlumes_smallP[n_1,],
                        dataLoc = d_1)  
)
int1_part1 = fun1(y = as.matrix(radPlumes_small@data[intersect(d_1, dS1_1),]),
                 data = radPlumes_smallP[intersect(d_1, dS1_1),],
                 newdata = radPlumes_smallP[intersect(n_1, nS1_1),])
int1_part2 = fun2(y = as.matrix(radPlumes_small@data[intersect(d_1, dS2_1),]),
                  data = radPlumes_smallP[intersect(d_1, dS2_1),],
                  newdata = radPlumes_smallP[intersect(n_1, nS2_1),])
int1_part3 = fun3(y = as.matrix(radPlumes_small@data[intersect(d_1, dS3_1),]),
                  data = radPlumes_smallP[intersect(d_1, dS3_1),],
                  newdata = radPlumes_smallP[intersect(n_1, nS3_1),])
int1_parts = matrix(nrow = length(n_1), ncol = ncol(int1))
int1_parts[is.element(n_1, nS1_1)] = int1_part1
int1_parts[is.element(n_1, nS2_1)] = int1_part2
int1_parts[is.element(n_1, nS3_1)] = int1_part3
expect_equivalent(int1, int1_parts)



# test inside of 'interpolate'
data(radioactivePlumes_local)
fun1 = idw0z 
fun2 = replaceDefault(krige0, newDefaults = list(formula = z ~ 1, beta = NA, model = vgm(10, "Sph", 50000, 5), ... = NA), type = "fun_interpolation.interpolate")[[1]]
fun3 = replaceDefault(krige0, newDefaults = list(formula = z ~ 1, beta = NA, model = vgm(100, "Sph", 500, 50), ... = NA), type = "fun_interpolation.interpolate")[[1]]
fun4 = replaceDefault(krige0, newDefaults = list(formula = z ~ 1, beta = NA, model = vgm(1000, "Sph", 50, 500), ... = NA), type = "fun_interpolation.interpolate")[[1]]

radioactivePlumes_local@locations@data$part = c(rep(1, 441), rep(2, 2001 - 441), rep(3, 4161 - 2001), rep(4, 7629 - 4161))
x = spplot(radioactivePlumes_local@locations, zcol = "part")
nS1 = list(1:441, 442:2001, 2002:4161, 4162:7629)
int1 = replaceDefault(fun = interpol0reg,
                         newDefaults = list(
                           newdataSplit = nS1,
                           dataSplit = 1:7629,
                           fun_interpolation = list(fun1, fun2, fun3, fun4)
                         ), type = "fun_interpolationSplit.interpolate")
loc1 = seq(1, 7629, 100)
interpolated1 = interpolate(simulations = radioactivePlumes_local,
                            locations = loc1,
                            values = 1,
                            fun_interpolation = int1[["fun"]])

radPl_local_points = as(radioactivePlumes_local@locations, "SpatialPointsDataFrame")    
names(radPl_local_points@data)[1] = "z"
radPl_local_values = matrix(getValues(subset(radioactivePlumes_local, kinds = 1)@values), byrow = TRUE, nrow = 7629)
interpolated1_part1 = fun1(y = radPl_local_values[loc1,],
                           data = radPl_local_points[loc1,],
                           newdata = radPl_local_points[nS1[[1]],])
interpolated1_part2 = fun2(y = radPl_local_values[loc1,],
                           data = radPl_local_points[loc1,],
                           newdata = radPl_local_points[nS1[[2]],])
interpolated1_part3 = fun3(y = radPl_local_values[loc1,],
                           data = radPl_local_points[loc1,],
                           newdata = radPl_local_points[nS1[[3]],])
interpolated1_part4 = fun4(y = radPl_local_values[loc1,],
                           data = radPl_local_points[loc1,],
                           newdata = radPl_local_points[nS1[[4]],])
interpolated1_parts = matrix(nrow = 7629, ncol = 876)
interpolated1_parts[nS1[[1]],] = interpolated1_part1 
interpolated1_parts[nS1[[2]],] = interpolated1_part2 
interpolated1_parts[nS1[[3]],] = interpolated1_part3 
interpolated1_parts[nS1[[4]],] = interpolated1_part4 

expect_equivalent(as.matrix(interpolated1), interpolated1_parts)
