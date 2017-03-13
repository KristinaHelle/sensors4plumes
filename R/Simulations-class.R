################################################################
#                    Simulations-class                        #
################################################################

setClass("Simulations",
         representation = representation(
           locations = "SpatialDataFrame",
           plumes = "data.frame",
           values = "Raster")
         )

setValidity("Simulations", function(object){
  if(is.null(object@locations)){ 
    return("'locations' missing")
  }else{
    nRow = nrow(object@values)
    nCol = ncol(object@values)
    if(!is.null(object@values)){
      nLocations = nrow(object@locations@data)
      if(nLocations != nRow){
        return(paste("The number of 'locations' (",nLocations,")
                     differs from the number of rows of 'values' (", nRow, ").
                     They have to agree, the rows of 'values' must belong to the locations that are represented by the rows of 'locations'."))
      }
    }
    if(!is.null(object@plumes)){
      nPlumes = nrow(object@plumes)  
      if(nPlumes != nCol){
        return(paste("The number of 'plumes' (", nPlumes, ")
                     differs from the number of columns of 'values' (", nCol, ").
                     They have to agree, the columns of 'values' must belong to the plumes that are represented by the rows of 'plumes'."))
      }
    }
    
    # object@values must have standard georeference to ensure that objects fit together, it has no meaning
    extent = c(object@values@extent@xmin, 
               object@values@extent@xmax, 
               object@values@extent@ymin, 
               object@values@extent@ymax)
    if (!all(extent == c(-90, 90, -90, 90))){
      return("Wrong extent of values, it has to be: -90, 90, -90, 90. To correct the extent of object@values of raster objects from files, change the '.grd' file.")
    }
    if (!identical(object@values@crs, CRS("+init=epsg:4326"))){
      return("Wrong projection of values, it has to be '+init=epsg:4326' (WGS 84). To correct the projection of object@values of raster objects from files, change the '.grd' file.")
    }
    
    return(TRUE)
  }
})

Simulations = function(locations,
                       plumes,
                       values){  
  new = new("Simulations",
            locations = locations,
            plumes = plumes,
            values = values)
  return(new)
}
