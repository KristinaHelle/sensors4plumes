##################################################################
#                subset.SpatialDataFrame                         #
##################################################################
# subsetSDF

subsetSDF = 
  function(x, locations, data = names(x@data), ...){}

subsetSDF.SpatialIndexDataFrame  = function(
  x, # SpatialIndexDataFrame
  locations, # integer; select the elements with this index value
  data = names(x@data), # integer/character; select these columns of data
  ...,
  grid # integer; select these elements of x@index
){
  if (!is.na(x@proj4string@projargs)){
    new_proj4string = x@proj4string
  }else{
    new_proj4string = CRS(as.character(NA))
  }
  # determine which index entries remain
  nI = length(x@index)
  keep_grid = rep(TRUE, nI)
  keep_index = rep(TRUE, nI)
  if(!missing(data)){
    new_data = x@data[,data, drop = FALSE]
  }else{
    new_data = x@data
  }
  if(!missing(grid)){
    keep_grid[!is.element(1:nI, grid)] = FALSE
  }
  if(!missing(locations)){
    keep_index[!is.element(x@index, locations)] = FALSE    
    if (any(table(locations) > 1)){
      warning("Some 'locations' have been selected multiple times, this is ignored.") 
    } 
    #if (!identical(sort(locations), locations)){
    #  warning("The order of 'locations' is ignored.")
    #}
  }
  new_index = x@index[keep_grid & keep_index]
  if (!missing(locations)){
    orderLocations = as.integer(order(unique(na.omit(intersect(locations, new_index)))))
    rankLocations = as.integer(rank(unique(na.omit(intersect(locations, new_index)))))
  }
  
  ## adapt index & data to changes from deleting and order of locations
  remaining_index = sort(unique(na.omit(new_index)))
  new_n = length(remaining_index)
  New_index = rep(NA, length(new_index))
  for(i in 1:new_n){
    if (missing(locations)){
      New_index[new_index == remaining_index[i]] = i         
    } else {
      New_index[new_index == remaining_index[i]] = orderLocations[i] 
    }
  }
  if (missing(locations)){
    new_data = new_data[remaining_index,, drop = FALSE]
  } else {
    new_data = new_data[remaining_index[rankLocations],, drop = FALSE]  
  }
  
  
  new = new("SpatialIndexDataFrame",
            index = New_index,
            data = new_data,
            bbox = matrix(0, nrow = 2, ncol = 2),
            proj4string = new_proj4string)
  return(new)
}
setMethod("subsetSDF", signature(x = "SpatialIndexDataFrame"), subsetSDF.SpatialIndexDataFrame)

subsetSDF.SpatialPointsDataFrame = 
  function(x, locations, data = names(x@data), ...){
    if (!missing(locations)){
      out = x[locations, data]
    } else {
      out = x[, data]
    }
    return(out)
  }
setMethod("subsetSDF", signature(x = "SpatialPointsDataFrame"), subsetSDF.SpatialPointsDataFrame)

subsetSDF.SpatialPixelsDataFrame = 
  function(x, locations, data = names(x@data), ...){
    if (!missing(locations)){
      out = x[locations, data]
    } else {
      out = x[, data]
    }
    return(out)
  }
setMethod("subsetSDF", signature(x = "SpatialPixelsDataFrame"), subsetSDF.SpatialPixelsDataFrame)

subsetSDF.SpatialPolygridDataFrame =  
  function (x, 
            locations,             # integer, index of index
            data = names(x@data),  # integer or character, index of data columns
            ...,
            coord_x, coord_y,      # numeric each (min, max) of grid x-/ y- coordinates (note: directions differ from grid_i, grid_j)
            grid_i, grid_j,        # integer, index of grid row / column
            grid_ij                # matrix same size as grid, logical
  ){
    # subset columns of data
    new_data = x@data[,data, drop = FALSE]
    
    if (all(missing(grid_i), missing(grid_j), missing(locations), missing(coord_x), missing(coord_y), missing(grid_ij))){
      new = x
      new@data = new_data 
    } else { # subset spatially
      indexNew = matrix(x@index, ncol = x@grid@cells.dim[1], nrow = x@grid@cells.dim[2], byrow = TRUE) 
      
      ## are subsetting parameters given?
      ##                              valid? 
      ##                                 result in empty grid?
      if(!missing(grid_i)){
        if(all(abs(grid_i - round(grid_i)) == 0)){
          if(max(grid_i) >= 1 & min(grid_i) <= x@grid@cells.dim[2]){
            selected = is.element(1:x@grid@cells.dim[2], grid_i)
            indexNew[!selected,] = NA  
            if(!all(is.element(grid_i, 1:x@grid@cells.dim[2]))){
              warning("grid_i contains indices beyond the number of rows of the grid, these are ignored.")
            }
          }else{# no overlap
            stop("grid_i are not row indices of the grid, the result is empty.")
          }
        }else{
          return("grid_i must be integer values; the given values are ignored for subsetting.")
        }
      }     
      
      if(!missing(grid_j)){
        if(all(abs(grid_j - round(grid_j)) == 0)){
          if(max(grid_j) >= 1 & min(grid_j) <= x@grid@cells.dim[1]){
            selected = is.element(1:x@grid@cells.dim[1], grid_j)
            indexNew[,!selected] = NA  
            if(!all(is.element(grid_j, 1:x@grid@cells.dim[1]))){
              warning("grid_j contains indices beyond the number of columns of the grid, these are ignored.")
            }
          }else{# no overlap
            stop("grid_j are not column indices of the grid, the result is empty.")
          }
        }else{
          return("grid_j must be integer values; the given values are ignored for subsetting.")
        }
      }     
      if(!missing(coord_x)){
        if(is.numeric(coord_x) & length(coord_x) == 2){
          if(max(coord_x) >= x@bbox[1,1] & min(coord_x) <= x@bbox[1,2]){
            coordIndex = c(
              floor((min(coord_x) - x@bbox[1,1])/x@grid@cellsize[1]) + 1,
              ceiling((max(coord_x) - x@bbox[1,1])/x@grid@cellsize[1]))         
            coordIndex[1] = max(coordIndex[1], 1)  
            coordIndex[2] = min(coordIndex[2], x@grid@cells.dim[1])
            selected = 1:x@grid@cells.dim[1] >= coordIndex[1] & 
              1:x@grid@cells.dim[1] <= coordIndex[2]
            indexNew[,!selected] = NA
            if(max(coord_x) > x@bbox[1,2] | min(coord_x) < x@bbox[1,1]){
              warning("coord_x exceeds the grid, the overlapping subset is used.")
            }
          }else{# no overlap
            stop("coord_x are not x-coordinats of the grid, the result is empty.")
          }      
        }else{
          warning("coord_x must be a vector of two numeric values; the given values are ignored for subsetting.")
        }
      }  
      if(!missing(coord_y)){
        if(is.numeric(coord_y) & length(coord_y) == 2){
          if(max(coord_y) >= x@bbox[2,1] & min(coord_y) <= x@bbox[2,2]){
            coordIndex = c(
              floor((min(coord_y) - x@bbox[2,1])/x@grid@cellsize[2]) + 1,
              ceiling((max(coord_y) - x@bbox[2,1])/x@grid@cellsize[2]))         
            coordIndex[1] = max(coordIndex[1], 1)  
            coordIndex[2] = min(coordIndex[2], x@grid@cells.dim[2])
            selected = 1:x@grid@cells.dim[2] >= coordIndex[1] & 1:x@grid@cells.dim[2] <= coordIndex[2]
            indexNew[rev(!selected),] = NA
            if(max(coord_y) > x@bbox[2,2] | min(coord_x) < x@bbox[2,1]){
              warning("coord_y exceeds the grid, the overlapping subset is used.")
            }
          }else{# no overlap
            stop("coord_y are not y-coordinats of the grid, the result is empty.")
          }      
        }else{
          warning("coord_y must be a vector of two numeric values; the given values are ignored for subsetting.")
        }
      } 
      if(!missing(grid_ij)){
        if(is.logical(grid_ij) && dim(t(grid_ij)) == x@grid@cells.dim){
          if (any(grid_ij)){
            indexNew[!grid_ij] = NA 
          } else {# empty grid
            stop("grid_ij contains only FALSE, the result is empty.")
          }
        }else{
          warning("grid_ij is invalid (not logical matrix of same size as grid), it is ignored for subsetting.")        
        }        
      }    
      
      if(!missing(locations)){
        if (any(table(locations) > 1)){
          warning("Some 'locations' have been selected multiple times, this is ignored.") 
        } 
#        if (!identical(sort(locations), locations)){
#          warning("The order of 'locations' is ignored.")
#        }
        if(is.numeric(locations)){
          if(any(is.element(locations, x@index))){
            indexGrid = matrix(is.element(x@index, locations),
                               nrow = x@grid@cells.dim[2],
                               ncol = x@grid@cells.dim[1],
                               byrow = TRUE)
            indexNew[!indexGrid] = NA
            if (!all(is.element(locations, x@index))){
              warning("'locations' contains values that do not occur, they are ignored")
            }
          }else{
            stop("The values given in 'locations' do not occur in the the index of x, the result is empty.")
          }
        }else{
          warning("The values in 'locations' are invalid (no integers), they are not used for subsetting.")
        }
      } 
      if (all(is.na(indexNew))){
        stop("Combining all subsetting parameters, the result is empty.")
      }
      
      # subsetting of the grid and adaption of index and data if necessary (i.e. if dropping rows or columns of the grid makes some values of index not occur anymore)  
      indexMargin1 = apply(X = is.na(indexNew), FUN = all, MARGIN = 1)
      indexMargin2 = apply(X = is.na(indexNew), FUN = all, MARGIN = 2)
      margin = rbind(c(min(which(indexMargin1 == FALSE)), max(which(indexMargin1 == FALSE))),
                     c(min(which(indexMargin2 == FALSE)), max(which(indexMargin2 == FALSE))))
      
#      if (margin[1,1] < margin[1,2] & margin[2,1] < margin[2,2]){
#        new_grid = x@grid[margin[1,1]:margin[1,2], margin[2,1]:margin[2,2]]

        new_grid = x@grid  
        new_grid@cellcentre.offset[1] = x@grid@cellcentre.offset[1] + 
                                        x@grid@cellsize[1] * (margin[2,1] - 1)
        new_grid@cells.dim[1] = as.integer(abs(diff(margin[2,])) + 1)
        new_grid@cellcentre.offset[2] = x@grid@cellcentre.offset[2] + 
                                        x@grid@cellsize[2] * (x@grid@cells.dim[2] - margin[1,2])
        new_grid@cells.dim[2] = as.integer(abs(diff(margin[1,])) + 1)
 
        new_spatialGrid = SpatialGrid(new_grid, x@proj4string)
        new_index = as.integer(t(indexNew[margin[1,1]:margin[1,2], margin[2,1]:margin[2,2]]))
        
        ## adapt index and data to changes from deleting
if (!missing(locations)){
  orderLocations = as.integer(order(unique(na.omit(intersect(locations, new_index)))))
  rankLocations = as.integer(rank(unique(na.omit(intersect(locations, new_index)))))
}
remaining_index = sort(unique(na.omit(new_index)))
new_n = length(remaining_index)
New_index = rep(NA, length(new_index))
for(i in 1:new_n){
  if (missing(locations)){
    New_index[new_index == remaining_index[i]] = i         
  } else {
    New_index[new_index == remaining_index[i]] = orderLocations[i] 
  }
}
if (missing(locations)){
  new_data = new_data[remaining_index,, drop = FALSE]  
} else {
  new_data = new_data[remaining_index[rankLocations],, drop = FALSE]
}

#        remaining_index = sort(unique(na.omit(new_index)))
#        new_n = length(remaining_index)
#        New_index = rep(NA, length(new_index))
#        for(i in 1:new_n){
#          New_index[new_index == remaining_index[i]] = i   
#        }
#        new_data = new_data[remaining_index,, drop = FALSE]
        new = new("SpatialPolygridDataFrame",
                  grid = new_grid,
                  index = New_index,
                  proj4string = x@proj4string,
                  bbox = bbox(new_spatialGrid),
                  data = new_data)     
#      } else { # result has only one row or column, must be transformed into SpatialPointsDataFrame
#        warning("The resulting grid would only have one row or column, therefore it is transformed to a SpatialPointsDataFrame.")
#        xPoints = as(x@grid, "SpatialPoints")
#        new_index = as.integer(t(indexNew))
#        xPoints = xPoints[!is.na(new_index)]
#        new_data = new_data[na.omit(new_index),, drop = FALSE] # copy data to the associated points
#        new = SpatialPointsDataFrame(coords = xPoints, data = new_data, proj4string = proj4string(x))
#      }  
    }
    return(new)
  }   
setMethod("subsetSDF", signature(x = "SpatialPolygridDataFrame"), subsetSDF.SpatialPolygridDataFrame)

subsetSDF.SpatialPolygonsDataFrame = 
  function(x, locations, data = names(x@data), ...){
    if (!missing(locations)){
      out = x[locations, data]
    } else {
      out = x[, data]
    }
    return(out)
  }
setMethod("subsetSDF", signature(x = "SpatialPolygonsDataFrame"), subsetSDF.SpatialPolygonsDataFrame)

subsetSDF.SpatialLinesDataFrame = 
  function(x, locations, data = names(x@data), ...){
    if (!missing(locations)){
      out = x[locations, data]
    } else {
      out = x[, data]
    }
    return(out)
  }
setMethod("subsetSDF", signature(x = "SpatialLinesDataFrame"), subsetSDF.SpatialLinesDataFrame)




