#####################################################
#            SpatialIndexDataFrame-class            #
#####################################################
setClass("SpatialIndexDataFrame",
         representation = representation(
           data = "data.frame",
           index = "integer",
           bbox = "matrix",
           proj4string = "CRS")
)

setValidity("SpatialIndexDataFrame", function(object){
  if(is.null(object@data)){
    return("'data' missing")
  }
  if(is.null(object@index)){
    return("'index' missing")
  }
  if( !is.null(object@data) & !is.null(object@index)){
    if(!all(sort(unique(object@index)) == seq(length = nrow(object@data)))){
      return("'index' does not consist of the numbers seq(length = nrow('data')).")
    }
  }  
  return(TRUE) 
})

SpatialIndexDataFrame = function(data, index
){
  if (!is.integer(index)){
    warning("'index' has the wrong class, it must be integer; it is ignored.")
  } else {
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
  }  
  new("SpatialIndexDataFrame",
      data = data,
      index = index,
      bbox = matrix(0, nrow = 2, ncol = 2),
      proj4string = CRS(as.character(NA)))
}
