############################################################
# test      testData.SpatialPolygridDataFrame              #
############################################################
# population_polygrid
# plumes_polygrid

data(population_polygrid)
spplot(population_polygrid, zcol = 5)

data(plumes_polygrid)
polygrid = cbind(population_polygrid, plumes_polygrid)
expect_equal(
  plumes_polygrid@data[,4],
  polygrid@data[,9]
)


