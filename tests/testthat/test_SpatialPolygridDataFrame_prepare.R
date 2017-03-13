###############################################################
# test SpatialPolygridDataFrame: prepare data                 #
###############################################################
index1 = as.integer(
  c( 6, 6, 7, 7, 8, 8,
     6, 6, 7, 7, 8, 8,
     5, 5, 1, 2, 9, 9,
     5, 5, 4, 3, 9, 9,
     12,12,11,11,10,10,
     12,12,11,11,10,10))
dataFrame1 = data.frame(a = 1:12 * 10, b = 2^(12:1))
dataFrame2 = data.frame(a = 12:1 * 10, b = 2^(1:12))
dataFrame3 = data.frame(a = 2 * 1:36, b = 0.5^{1:36}, c = runif(36))

grid1 = GridTopology(c(1, 1), c(2, 2), c(6, 6))
grid2 = GridTopology(c(1, 1), c(0.99, 0.43), c(6, 14))
grid3 = GridTopology(c(13, 13), c(2, 2), c(6, 6))
grid4 = GridTopology(c(1, 1), c(1, 4), c(6, 6))
spatialGrid1 = SpatialGrid(grid1, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
spatialGrid2 = SpatialGrid(grid1, CRS("+init=epsg:4267"))
spatialPoints1 = spatialGrid1
gridded(spatialPoints1) = FALSE
coordinates1 = spatialPoints1@coords
coordinates2 = coordinates1 + runif(length(coordinates1))
coordinates3 = coordinates(grid4)
data(plumes_polygrid)
coordinates4 = coordinates(plumes_polygrid)
spatialGridDataFrame1 = SpatialGridDataFrame(grid = spatialGrid1, 
                                             data = dataFrame1[index1,])
spatialPointsDataFrame1 = SpatialPointsDataFrame(coords = spatialPoints1, 
                                                 data = dataFrame1[index1,])
spatialPointsDataFrame2 = SpatialPointsDataFrame(
  coords = coordinates2,
  data = dataFrame1[index1,]
)
point1 = data.frame(x = 0, y = 0)
