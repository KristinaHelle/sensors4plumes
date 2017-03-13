######################################################
#           SpatialPolygridDataFrame-class           #
######################################################

setClass("SpatialPolygridDataFrame",
         representation = representation(
           data = "data.frame",
           index = "integer",
           grid = "GridTopology",
           bbox = "matrix",
           proj4string = "CRS"))


setValidity("SpatialPolygridDataFrame", function(object){
  if(is.null(object@data)){
    return("'data' missing")
  }
  if(is.null(object@index)){
    return("'index' missing")
  }
  if(is.null(object@grid)){ 
    return("'grid' missing")
  }
  if(is.null(object@proj4string)){
    return("'proj4string' missing")
  }

  if(is.null(object@bbox)){
    return("'bbox' missing")
  } else {  
    # check bbox by spatial grid
    spatialGrid = SpatialGrid(grid = object@grid, proj4string = object@proj4string)
    if (! identical(object@bbox, bbox(spatialGrid))){
      warning(paste0("'bbox' is not correct, it must be the correct bbox of 'grid':",
            bbox(spatialGrid)))
    }
  }
  
  if(!is.null(object@grid@cells.dim) & !is.null(object@data) & !is.null(object@index)){
    cellsDim = object@grid@cells.dim
    if(prod(cellsDim) != length(object@index)){
      return(paste0("Number of cells in the grid (",prod(cellsDim),") 
                    does not fit length of the index (", length(object@index), ")."))
    }
    
    indexUnique = sort(unique(object@index))
    textIndex = paste0("'index' must consist of the rownumbers of 'data' 1:", nrow(object@data), ".")
    if (length(indexUnique) == nrow(object@data)){
      if (!all(indexUnique == seq(length = nrow(object@data)))){
        return(textIndex)      
      }
    } else {
      return(textIndex)
    }
  }
})

# function to generate objects of class SpatialPolygridDataFrame
# based on a spatial grid or on a grid topology (and a projection)
SpatialPolygridDataFrame = function(
  data, 
  index = 1:nrow(data),
  grid, 
  proj4string
){
  if(is(grid, "SpatialGrid")){
    proj4string_new = CRS(sp::proj4string(grid)) 
    grid = grid@grid
  } else {
    proj4string_new = CRS(as.character(NA))
  }
  if (!missing(proj4string)){
    proj4string_new = proj4string
  }
  
  if (!missing(index)){
    if (!missing(grid)){
      if (prod(grid@cells.dim) != length(index)){
        stop("length of 'index' does not fit number of cells in 'grid'.")
      }
    }
    if (!is.integer(index)){
      stop("'index' must be of class integer.")
    }
    indexNeg = index <= 0 
    if (any(indexNeg, na.rm = TRUE)){
      warning("'index' contains negative values, they are replaced by 'NA'.")
      index[indexNeg] = NA
    }
    if (!missing(data)){
      indexTooHigh = index > nrow(data)
      if (any(indexTooHigh, na.rm = TRUE)){
        warning("'index' contains values that exceed the number of 'data', they are replaced by 'NA'.")
        index[indexTooHigh] = NA
      } 
    }  
  }
  # if not all possible values occur in index, adjust it, some data may be lost
  remaining_index = sort(unique(na.omit(index)))
  new_n = length(remaining_index)
  if(new_n < max(index, na.rm = TRUE)){
    new_index = rep(NA, length(index))
    for(i in 1:new_n){
      new_index[index == remaining_index[i]] = i   
    }
    index = new_index
  }
  if (!missing(data)){
    if (new_n < nrow(data)){
      warning("Some of the 'data' is ignored as it is not associated to 'index' values.")
    }  
    data = data[remaining_index,, drop = FALSE]     
  } else {
    data = data.frame(Index = 1:new_n)  
  }
  

  
  # generate spatial grid
  spatialGrid = SpatialGrid(grid = grid, proj4string = proj4string_new)
  bbox_new = bbox(spatialGrid)

  
  new("SpatialPolygridDataFrame", 
      data = data,
      index = index,
      grid = grid,
      bbox = bbox_new,
      proj4string = proj4string_new)
}
