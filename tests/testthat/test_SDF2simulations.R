########################################################################
# test          SDF2Simulations                                        #
########################################################################

# small example, no 'indices' given (generate from Simulations)
data(population_polygrid)  
population_simulations = SDF2simulations(population_polygrid)
population_polygrid2 = extractSpatialDataFrame(population_simulations)  
expect_equal(
  dename(dename(population_polygrid2@data), kind = "rownames"),
  dename(dename(population_polygrid@data), kind = "rownames"))

# small example with multiple-row 'indices'
Indices1 = matrix(c(1,3,5, 6,8,10), byrow = TRUE, nrow = 2, dimnames = list(c("kind1", "kind2"), 1:3))
sSPixelsDF1 = SDF2simulations(eSPixelsDF1, indices = Indices1)
expect_equal(
  sSPixelsDF1@locations@grid, eSPixelsDF1@grid)
expect_equal(
  sSPixelsDF1@locations@proj4string, eSPixelsDF1@proj4string)
expect_equal(
  sSPixelsDF1@locations@data,
  eSPixelsDF1@data[,setdiff(1:10, as.integer(Indices1))])
expect_equal(
  getValues(sSPixelsDF1@values)[1:3,1], 
  c(eSPixelsDF1@data[1,Indices1[1,1]],
    eSPixelsDF1@data[1,Indices1[1,2]],
    eSPixelsDF1@data[1,Indices1[1,3]]))
expect_equal(
  names(sSPixelsDF1@values),
  dimnames(Indices1)[[1]])


# big example with multiple-row 'indices'
sRadioactivePlumes_area_725_maxdose = SDF2simulations(RadioactivePlumes_area_725_maxdose)
expect_is(sRadioactivePlumes_area_725_maxdose, "Simulations")

# error: missing|wrong values in 'indices' 
Indices2 = matrix(c(1,3,NA, 6,8,10), byrow = TRUE, nrow = 2, dimnames = list(c("kind1", "kind2"), 1:3))
expect_error(SDF2simulations(eSPixelsDF1, indices = Indices2))
Indices3 = matrix(c(1,3,5, 6,8,12), byrow = TRUE, nrow = 2, dimnames = list(c("kind1", "kind2"), 1:3))
expect_error(SDF2simulations(eSPixelsDF1, indices = Indices3))

