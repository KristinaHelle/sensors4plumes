#######################################################################
# test        spplot.SpatialPolyridDataFrame                          #
#######################################################################
spplot(SPolygridDFB2, zcol = "a", col.regions = rainbow(20), 
       sp.layout = list("sp.points", coordinates(SPolygridDFB2), col = 1))

expect_error(# 1 spatial point, cannot be transformed into SpatialPolygridDataFrame
  spplot(points2polygrid(point1))  
)

