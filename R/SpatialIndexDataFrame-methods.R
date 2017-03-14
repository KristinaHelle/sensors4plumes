###############################################################
#                 SpatialIndexDataFrame-methods               #
###############################################################
# coerce: data.frame ->

# coordinates
# proj4string
# spplot
# cbind
# rbind
# subsetSDF => subsetSDF
# areaSDF => areaSDF

setAs("data.frame", "SpatialIndexDataFrame",
      function(from){
        to = new("SpatialIndexDataFrame",                 
                 data = from, 
                 index = seq(length = nrow(from)),
                 bbox = matrix(0, nrow = 2, ncol = 2),
                 proj4string = CRS(as.character(NA)))   
        return(to)
      }
)


setMethod("coordinates", signature(obj = "SpatialIndexDataFrame"),
          function(obj){ 
            warning("SpatialIndexDataFrame objects do not have coordinates.")
            out = NULL
            return(out)
          }
)

setMethod("proj4string", signature(obj = "SpatialIndexDataFrame"),
          function(obj){ 
            warning("SpatialIndexDataFrame objects do not have a projection.")
            out = as.character(NA)
            return(out)
          }
)

setMethod("bbox", signature(obj = "SpatialIndexDataFrame"),
          function(obj){ 
            warning("SpatialIndexDataFrame objects have no spatial extent.")
            bbox = obj@bbox
            return(bbox)
          }
)

setMethod("spplot", signature(obj = "SpatialIndexDataFrame"),
          function(obj, ...){ 
            values = obj@data[obj@index,]
            values$x = 1
            values$y = rev(seq(length = length(values[[1]])))  
            #values$index = obj@index
            gridded(values) = ~ x + y
            #values@grid@cellsize["x"] = values@grid@cells.dim["y"]/10
            print(spplot(values, ...))
          }  
)

setMethod("length", signature(x = "SpatialIndexDataFrame"),
          function(x){ 
            length = nrow(x@data)
          }  
)

cbind.SpatialIndexDataFrame = function(...){ 
  dots = list(...)
  n = length(dots)
  if (n >= 1){
    newObj = dots[[1]]   
  }else{
    stop("No object(s) given! Returns NULL.")
    newObj = NULL
  }
  if (n >= 2){
    for (l in 2:n){
      if(!identical(newObj@index, dots[[l]]@index)){
        stop("SpatialIndexDataFrames cannot be combined as their 'index' differ.")        
      } else {
        newObj = new("SpatialIndexDataFrame",
                     index = newObj@index,
                     proj4string = newObj@proj4string,
                     bbox = newObj@bbox,
                     data = cbind(newObj@data, dots[[l]]@data)
        )
        
      }
    }
  }
  return(newObj)
}  
   
# rbindSpatialIndexDataFrame = function(...){ 
#   dots = list(...)
#   n = length(dots)
#   if (n >= 1){
#     newObj_data = dots[[1]]@data
#     newObj_index = dots[[1]]@index
#   }else{
#     stop("No object(s) given! Returns NULL.")
#     newObj = NULL
#   }
#   if (n >= 2){
#     for (l in 2:n){
#       n_data = nrow(newObj_data)
#       newObj_data = rbind(newObj_data, dots[[l]]@data)
#       newObj_index = c(newObj_index, dots[[l]]@index + n_data)
#     }
#   }  
#   newObj = new("SpatialIndexDataFrame",
#                      index = newObj_index,
#                      proj4string = dots[[1]]@proj4string,
#                      bbox = dots[[1]]@bbox,
#                      data = newObj_data)
#   return(newObj)
# }  




