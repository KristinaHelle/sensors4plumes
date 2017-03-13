####################################################
#                 SpatialDataFrame-class           #
####################################################

setClassUnion(name = "SpatialDataFrame", 
              members = c("SpatialIndexDataFrame", 
                          "SpatialPointsDataFrame", 
                          "SpatialPixelsDataFrame", 
                          #"SpatialGridDataFrame", 
                          "SpatialPolygridDataFrame", 
                          "SpatialPolygonsDataFrame", 
                          "SpatialLinesDataFrame"))


is.SpatialDataFrame = function(SpatialDataFrame){
  classSDF = class(SpatialDataFrame)
  if (is.element(classSDF, c("SpatialIndexDataFrame", 
                             "SpatialPointsDataFrame", 
                             "SpatialPixelsDataFrame", 
                             "SpatialPolygridDataFrame", 
                             "SpatialPolygonsDataFrame", 
                             "SpatialLinesDataFrame"))){
    out = TRUE
  } else {
    out = FALSE
  }
  return(out)
}








