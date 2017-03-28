##################################################################
#              SpatialPolygridDataFrame-methods                  #
##################################################################
# these methods are based on these functions:
# points2polygrid
# polygrid2grid

# points2polygrid-----------------------------------------------
points2polygrid = function(
  points,     # original coordinates: data.frame, or matrix with coordinates in columns 1 and 2; or SpatialPoints
  grid,  # target grid: "GridTopology" or "SpatialGrid" (sp)
  index, # target index: integer of length = ncells in grid, values from 1:nrow(points)
  data,  # data to be assigned, each row belongs to the coordinates in the same row of points
  tolerance = signif(mean(diff(t(bbox(points))))/1000, digits = 1) # point coordinates are rounded to this dist, smaller differences are ignored
  #tolerance = sqrt(.Machine$double.eps)
){
#  require(sp)
# message(paste0("tolerance is ", tolerance))
  parameters = c(!missing(points), !missing(grid), !missing(index), !missing(data))

  # proj4string
  # test if parameters have the correct form, else ignore
  if (all(parameters[1:2])){
    if (is(grid, "SpatialGrid") & is(points, "SpatialPoints")){
      if (proj4string(grid) != proj4string(points)){
        warning(paste0("The projections of 'grid' and 'points' differ, the one of the points (",
                       proj4string(points), ") is used."))
      }
    }
  }
  proj4string_new = as.character(NA)

  if (parameters[2]){
    if (is(grid, "SpatialGrid")){
      proj4string_new = proj4string(grid)
      grid = grid@grid
    }
    if (!is(grid, "GridTopology")){
      warning("'grid' has the wrong class, it must be SpatialGrid or GridTopology; it is ignored.")
      parameters[2] = FALSE
    }
  }

  if (parameters[1]){
    if (is(points, "SpatialPoints")){
      proj4string_new = proj4string(points)
      points = coordinates(points)
    } else {
      if(!is.element(class(points), c("matrix", "data.frame"))){
        warning("'points' has the wrong class, it must be SpatialPoints, matrix, or data.frame; it is ignored.")
        parameters[1] = FALSE
      } else {
        if(dim(points)[2] < 2){
          warning("'points' must have at least two columns, it is ignored.")
          parameters[1] = FALSE
        } else {
          if (!is.numeric(points[[1]]) | !is.numeric(points[[2]])){
            warning("The columns of 'points' must be numeric; 'points' is ignored.")
            parameters[1] = FALSE
          } else {
            points = points[,1:2]
          }
        }
      }
    }
  }

  if (all(parameters[2:3])){
    if(prod(grid@cells.dim) != length(index)){
      warning(paste("Number of cells in the grid (",prod(grid@cells.dim),
                    ") does not fit length of the index (", length(index),
                    "), 'index' is ignored."))
      parameters[3] = FALSE
    }
  }
  if (parameters[3]){
    if (!is.integer(index)){
      warning("'index' has the wrong class, it must be integer; it is ignored.")
      parameters[3] = FALSE
    } else {
      indexNeg = index <= 0
      if (any(indexNeg, na.rm = TRUE)){
        warning("'index' contains negative values, they are replaced by 'NA'.")
        index[indexNeg] = NA
      }
    }
  }

  # decide, if parameters are sufficient and which ones to use
  if (all(parameters[1:3] == c(FALSE, FALSE, FALSE)) |
        all(parameters[1:3] == c(FALSE, FALSE, TRUE)) |
        all(parameters[1:3] == c(FALSE, TRUE, FALSE))){
    stop("Not enough parameters given, either 'points' or 'index' and 'grid' needed.")
  }
  if(all(parameters[1:3])){
    warning("As 'index' and 'grid' are given, 'points' is ignored.")
    parameters[1] = FALSE
  }

  # if no grid given, determine it from points
  if(!parameters[2]){

    # occuring coordinates
    orig_x = sort(unique(points[,1]))
    orig_y = sort(unique(points[,2]))
    x = orig_x#sort(unique(orig_x))
    y = orig_y#sort(unique(orig_y))

    hasExtent = c(length(x) > 1, length(y) > 1)
    if (all(!hasExtent)){
      stop("There is only one 'point', it is impossible to determine the 'grid'.")
    }
    if (sum(hasExtent) == 1){
      warning("All 'points' lie on one line, grid resolution is determined from the other coordinate.")
    }

    # test points for regularity, if regular use dist between points as cellsize
    cellsize_x = NA
    cellsize_y = NA
    if (hasExtent[1]){
      uniq_dx = diff(x)
      min_dx = min(uniq_dx)
      if(max(uniq_dx%%min(uniq_dx)) == 0){
        cellsize_x = min_dx
#        message("cellsize_x from original")
        if (!hasExtent[2]){
          cellsize_y = cellsize_x
#          message("cellsize_y from x original")
        }
      }
    }
    if (hasExtent[2]){
      uniq_dy = diff(y)
      min_dy = min(uniq_dy)
      if(max(uniq_dy%%min(uniq_dy)) == 0){
        cellsize_y = min_dy
#        message("cellsize_y from original")
        if (!hasExtent[1]){
          cellsize_x = cellsize_y
#          message("cellsize_x from y original")
        }
      }
    }
    if (all(!is.na(c(cellsize_x, cellsize_y)))){
#      message("regular")
    } else {# not regular
      # round to toleranceAbs
      if (hasExtent[1] & is.na(cellsize_x)){
        x = unique(round(orig_x/tolerance) * tolerance)
      }
      if (hasExtent[2] & is.na(cellsize_y)){
        y = unique(round(orig_y/tolerance) * tolerance)
      }
      hasExtent = c(length(x) > 1, length(y) > 1)
      if (sum(hasExtent) == 1){
        warning("After rounding to tolerance all 'points' lie on one line, grid resolution is determined from the other coordinate.")
      }
      if (all(!hasExtent)){
        stop("When rounded to tolerance there is only one 'point', it is impossible to determine the grid.")
      }
      # test for regularity after rounding
      if (hasExtent[1] & is.na(cellsize_x)){
        uniq_dx = diff(x)
        min_dx = min(uniq_dx)
        if(max(uniq_dx%%min(uniq_dx)) == 0){
          cellsize_x = min_dx
#          message("cellsize_x from rounded")
          if (!hasExtent[2]){
            cellsize_y = cellsize_x
#            message("cellsize_y from rounded x")
          }
        }
      }
      if (hasExtent[2] & is.na(cellsize_y)){
        uniq_dy = diff(y)
        min_dy = min(uniq_dy)
        if(max(uniq_dy%%min(uniq_dy)) == 0){
          cellsize_y = min_dy
#          message("cellsize_y from rounded")
          if (!hasExtent[1]){
            cellsize_x = cellsize_y
#            message("cellsize_x from rounded y")
          }
        }
      }
    }
    if (all(!is.na(c(cellsize_x, cellsize_y)))){
#      message("regular after rounding to tolerance")
    } else {
      # determine grid resolution
      if (hasExtent[1] & is.na(cellsize_x)){
        # distances between coordinates, as multiples of 'tolerance'
        dx = as.integer(sort(unique(round(diff(x)/tolerance))))
        ## prime factors of distances as multiples of 'tolerance'
        primx = factorize(dx)
        # common divisor
        gcd_x = 1
        factor_x = unique(primx[[1]])
        if (! is.element(1, unlist(primx))){    # 1 (occurs only if some dx are 1)
          for (i in seq(along = factor_x)){
            j = 1
            while(all(dx%%factor_x[i]^j == 0)){
              gcd_x = gcd_x * factor_x[i]
              j = j + 1
            }
          }
        }
        # cellsize
        cellsize_x = gcd_x * tolerance
#        message("cellsize_x from common divisor")
        if (!hasExtent[2]){
          cellsize_y = cellsize_x
#          message("cellsize_y from common divisor of x")
        }
      }

      if (hasExtent[2] & is.na(cellsize_y)){
        dy = as.integer(sort(unique(round(diff(y)/tolerance))))
        primy = factorize(dy)
        gcd_y = 1
        factor_y = unique(primy[[1]])
        if (! is.element(1, unlist(primy))){
          for (i in seq(along = factor_y)){
            j = 1
            while(all(dy%%factor_y[i]^j == 0)){
              gcd_y = gcd_y * factor_y[i]
              j = j + 1
            }
          }
        }
        cellsize_y = gcd_y * tolerance
#        message("cellsize_y from common divisor")
        if (!hasExtent[2]){
          cellsize_y = cellsize_x
#          message("cellsize_x from common divisor of y")
        }
      }
    }

    if (hasExtent[1]){
      # grid dimension (number of cells)
      xn = as.integer(round(diff(range(x))/cellsize_x + 1))
    } else {
      cellsize_x = cellsize_y
      xn = 1
    }

    if (hasExtent[2]){
      yn = as.integer(round(diff(range(y))/cellsize_y + 1))
    } else {
      cellsize_y = cellsize_x
      yn = 1
    }

    # grid
    grid = GridTopology(c(x[1], y[1]),
                        c(cellsize_x, cellsize_y), c(xn, yn))
  }

  # if no index given, determine it as nearest neighbours
  if(!parameters[3]){
  #    require(FNN)
    gridPoints = coordinates(grid)
    neighbours = get.knnx(points, gridPoints, k = 1, algorithm = "kd_tree")
    index = neighbours$nn.index[,1]
    if (!all(is.element(1:nrow(points), index))){
      warning("Some of the points were merged, data$Index gives the indices of the points that were used.")
    }
  }

  # check, if data fits points / index
  if (parameters[4]){
    indexTooHigh = index > nrow(data)
    if (any(indexTooHigh, na.rm = TRUE)){
      if (parameters[3]){
        warning("'index' contains values that exceed the number of 'data', they are replaced by 'NA'.")
      }
      if (parameters[1]){
        warning(paste0("There are more 'points' (", max(index),
                       ") than rows of 'data' (", nrow(data),
                       "). Data is assigned to the first points, 'NA' to the extra points."))
      }
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
  if (parameters[4]){
    if (new_n < nrow(data)){
      warning("Some of the 'data' is ignored as it is not associated to points/index values.")
    }
    data = data[remaining_index,, drop = FALSE]
  } else {
    data = data.frame(Index = 1:new_n)
  }
#  } else {
#    if (parameters[4] == TRUE){
#      data = data
#    } else {
#      data = data.frame(Index = 1:max(index))
#    }
#  }

  out = SpatialPolygridDataFrame(
    grid = SpatialGrid(grid, CRS(proj4string_new)),
    index = index,
    data = data)
  return(out)
}

# polygrid2grid------------------------------------------------
polygrid2grid = function(
  obj,
  zcol = NA,  # ncol(obj@data)
  returnSGDF = TRUE,
  geoTiffPath
  ){
  if (all(is.na(zcol))){
    zcol = names(obj@data)
  }
  if (!is.character(zcol)){
    zcol = names(obj@data)[zcol]
  }
  zcol = c(zcol, "index")
  cells = SpatialGrid(obj@grid, obj@proj4string)
  fullgrid(cells) = FALSE
  cells = as(cells, "data.frame")
#  names(cells) = c("X", "Y")
  cells$index = obj@index
  imageValue = obj@data
  imageValue$index = seq(along = imageValue[[1]])
  # NA in index seems to cause cells to be deleted, replace by -1 (did not change the warning)
#  cells$index[is.na(cells$index)] = -1
#  imageValue = rbind(imageValue, imageValue[1,])
#  imageValue[nrow(imageValue),] = NA
#  imageValue[nrow(imageValue), "index"] = -1
  imageValueAll = merge(x = cells, y = imageValue, by = "index", all = TRUE, sort = FALSE)
#  imageValueAll$index[imageValueAll$index == -1] = NA
#  imageValueAll$index = NULL # changed!
  coordinates(imageValueAll) = 2:3 #~ x + y # (x,y does not work if imageValue has columns of this name)
  proj4string(imageValueAll) = obj@proj4string
  if (nrow(imageValueAll) <= 1){
    warning("Cannot derive grid parameters from a single point, SpatialPointsDataFrame returned, no geoTiff generated!")
    out = imageValueAll[,zcol, drop = FALSE]
  } else {
    gridded(imageValueAll) = TRUE
    fullgrid(imageValueAll) = TRUE
#    imageValueAll@data = imageValueAll@data[,zcol, drop = FALSE]
#    out = imageValueAll
    out = imageValueAll[,,zcol]

    if(!missing(geoTiffPath)){
      writeGDAL(dataset = imageValueAll[,,zcol], drivername = "GTiff", fname = paste(geoTiffPath, ".tif", sep = ""))
    }
  }

  if(returnSGDF){
    return(out)
  }else{
    return(TRUE)
  }
}


