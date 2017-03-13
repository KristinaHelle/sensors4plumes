######################################################################
# test spplotLog                                                     #
######################################################################
# test is to be done "by hand" as results need visual checking

data(SIndexDF)
data(SPointsDF)
data(SPixelsDF)
data(SPolygonsDF)
data(SLinesDF)
data(SPolygridDF)

# relative difference small => base = 1
SIndexDF
# relative difference medium -> base = 2
SIndexDF@data$d = c(0.01, 69, 288)
# relative difference big -> base = 10
SIndexDF@data$e = c(0.0000091, 46766, 286452348)
# zero, negative
SIndexDF@data$f = c(0, -1, 1)

# several or 1 attribute
spplotLog(SIndexDF, zcol = c("a", "b")) # base = 1
spplotLog(SIndexDF, zcol = "a") # base = 1

spplotLog(SIndexDF, zcol = c("a", "b", "f")) # base = 1, include <0, 0
spplotLog(SIndexDF, zcol = c("a", "b", "d")) # base = 2
spplotLog(SIndexDF, zcol = c("b", "e")) # base = 10, tickDiff = 2

# replace0
spplotLog(SIndexDF, zcol = c("a", "b", "c", "f"), replace0 = TRUE) # base = 1
spplotLog(SIndexDF, zcol = c("b", "e", "f"), replace0 = TRUE) # base = 10, tickDiff = 2
spplotLog(SIndexDF, zcol = c("a", "b", "d", "f"), replace0 = TRUE) # base = 2

# base
spplotLog(SIndexDF, zcol = c("a", "b", "d", "e", "f"), base = 5) # tickDiff = 2
spplotLog(SIndexDF, zcol = c("a", "b", "d", "e", "f"), replace0 = TRUE, base = 5) # tickDiff = 2
spplotLog(SIndexDF, zcol = c("a", "b", "d", "e", "f"), replace0 = TRUE, base = 50) # tickDiff = 1
spplotLog(SIndexDF, zcol = c("a", "b"), base = 5) # no nice ticks

# minNonzero
spplotLog(SIndexDF, zcol = c("a", "b", "d", "e", "f")) # tickDiff = 2
spplotLog(SIndexDF, zcol = c("a", "b", "d", "e", "f"), minNonzero = 1e-5) # tickDiff = 1, < 1e-5 as NA
spplotLog(SIndexDF, zcol = c("a", "b", "d", "e", "f"), minNonzero = 1e-5, replace0 = TRUE) #

# nTicks
spplotLog(SIndexDF, zcol = c("a", "b", "d", "e", "f"), base = 5, nTicks = 20) # tickDiff = 1
spplotLog(SIndexDF, zcol = c("a", "b", "d", "e", "f"), base = 5, nTicks = 50) # tickDiff < 1

# forward parameters to spplot
spplotLog(SIndexDF, zcol = c("a", "b", "d", "e", "f"), col.regions = rainbow, main = "test") # tickDiff < 1

# other objects
SPixelsDF@data$a = c(-5, 0, 0.003, 10, 57, 320, 444, 1000, 10000)
spplotLog(SPixelsDF, replace0 = TRUE, minNonzero = 0.05, col.regions = grey.colors,
          sp.layout = list("sp.points", SPointsDF[1,], col = 3, pch = 1, cex = 3))

spplotLog(SPolygonsDF, minNonzero = 15, replace0 = TRUE) 
spplotLog(SLinesDF, minNonzero = 15, replace0 = TRUE)
spplotLog(SPolygridDF, , minNonzero = 7, replace0 = TRUE, base = 2, col.regions = terrain.colors)
# SpatialGridDataFrame?

# SPointsDF
SPointsDF@data$a = c(0.01, 0.05, 1, 5, 100, 500, 0, 10000, 50000)
spplotLog(SPointsDF) # 0 as NA
spplotLog(SPointsDF, replace0 = TRUE) 
spplotLog(SPointsDF, replace0 = TRUE, minNonzero = 0.1) 
spplotLog(SPointsDF, replace0 = TRUE, minNonzero = 0.1, zcol = "a", col.regions = cm.colors(4), 
          sp.layout = list("sp.points", SPointsDF[1,], col = 3, pch = 1, cex = 3)) 

