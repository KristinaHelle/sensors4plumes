##################################################################
#              areaSDF.SpatialDataFrame                          #
##################################################################
# areaSDF

areaSDF = function(x){}

setMethod("areaSDF", signature(x = "SpatialIndexDataFrame"), function(x){
   warning(paste0("Objects of class ", class(x), " do not have area, 0 returned."))
   area = rep(0, nrow(x@data))
   return(area)
   }
)

setMethod("areaSDF", signature(x = "SpatialPointsDataFrame"), function(x){
  warning(paste0("Objects of class ", class(x), " do not have area, 0 returned."))
  area = rep(0, nrow(x@data))
  return(area)
  }
)

setMethod("areaSDF", signature(x = "SpatialPixelsDataFrame"), function(x){
  area = rep(prod(x@grid@cellsize), nrow(x@data)) 
  return(area)
  }
)

setMethod("areaSDF", signature(x = "SpatialPolygridDataFrame"), function(x){ 
  area = as.numeric(table(x@index) * prod(x@grid@cellsize)) 
  return(area)
  }
)

setMethod("areaSDF", signature(x = "SpatialPolygonsDataFrame"), function(x){ 
  area = sapply(X = x@polygons, FUN = slot, "area") 
  return(area)
  }
)

setMethod("areaSDF", signature(x = "SpatialLinesDataFrame"), function(x){ 
  area = SpatialLinesLengths(x)
  return(area)
  }
)
