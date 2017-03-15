# interpolation of all plumes of a Simulations object, 
## based on values in some common locations 
## solve: data not in memory

## computation of errors: extra functions, use via simulationsApply

## assumption: interpolation function: works with data, newdata, y as krige0, idw0
## multiplity and order of locations are ignored
## allow multiple plumes? no, if change of plumes is desired, subsample first 
###  (is needed later anyway for comparison); 
# similar to values
# -> simulations input must be the one we want interpolations for, no subsetting inside
## output always as raster


## computation plume (block)-wise, writing to raster location-wise?
### solution1: first save transposed (row ~ plume), then transpose back
# test
#y = stack("bigSimulations_finaldose.grd",
#          "bigSimulations_maxdose.grd",
#          "bigSimulations_time.grd")
#z = t(y) # works, but slow
### solution2: save immediately in final form:
#### leave unknown part of rows empty
#### when continuing to write to a row, first read, then combine and write

# parameters of errorInterpolation  
# simulations, # "Simulations" 
# locations, # "integer", indices of simulations@locations to be used as input for interpolation (multiples ignored)
# plumes, # "integer", indices of simulations@plumes to be used 
# values, # 'integer' or 'character'; layers to be kept
# #  keepData = FALSE, # if all values at these locations and plumes are to be kept (may not work if too much data)
# interpol, # interpolation function: krige0(), idw0(), or a function working the same way
# error_value, # (input, output, initial) must be cumulative (it can be computed blockwise, using the result as the new initial)
# error_map, # (input, output)
# tmpfile_error = "tmp_error", # filename, if results don't fit into memory: tmp_simApply_global.grd etc.; if FALSE: nothing saved to file, stop if it would be necessary
# tmpfile_interpol = FALSE, # filename (if does not fit into memory), FALSE ~ interpolation not saved
# overwrite = TRUE,
# #  na.rm = FALSE,
# ... # 

interpolate = function(
  simulations,
  locations,
  kinds = 1,
  fun_interpolation,
  tmpfile = "tmp_interpolate",
  overwrite = FALSE,
  chunksize = 1e+7
){
  # check class of simulations@locations, if possible turn into Spatial[sp] object
  if (is(simulations@locations,"SpatialIndexDataFrame")){
    stop("No interpolation possible as 'simulations' have 'locations' of class 'SpatialIndexDataFrame'.")
  }
  if (is(simulations@locations,"SpatialPolygridDataFrame")){
    newData = as(simulations@locations, "SpatialPointsDataFrame")
  } else {
    newData = simulations@locations
  } 
  names(newData)[1] = "z"
  # keep only first layer of simulations
  if (nlayers(simulations@values) > 1){
    simulations@values = subset(simulations@values, kinds[1], drop = FALSE)
  }
  
  # check function
  ## function with regions split?
  if(all(is.element(c("dataLoc", "newdataLoc"), names(formals(fun_interpolation))))){
    regFun = TRUE
  } else {
    regFun = FALSE
  }
  
  if (!regFun){
    funTested = replaceDefault(fun_interpolation, 
                               newDefaults = list(newdata = newData),
                               type = "fun_interpolation.interpolate")    
  } else {
    funTested = replaceDefault(fun_interpolation, 
                               newDefaults = list(newdata = newData),
                               type = "fun_interpolationSplit.interpolate")
  }

  fun_interpol = funTested[["fun"]]
  if (!funTested[["accept"]]){
    stop("No valid 'fun_interpolation'.")
  }
  # clean 'locations'
  nL = nLocations(simulations)
  if (!missing(locations)){
    locationsIn = locations > 0 & locations <= nL
    locations[!locationsIn] = NA
    locations = locations[!is.na(locations)]
    locations = sort(unique(locations))
    if (length(locations) < length(locationsIn)){
      warning("Some 'locations' are ignored because they are multiple, out of range or NA.")
    }
  }
  
  # check, what fits into memory
  #  bs_subset = blockSize(raster(nrow = length(locations), ncol = nPlumes(simulations)), 
  #                        n = nlayers(simulations@values), 
  #                        minblocks = 1, chunksize = chunksize)
  bs_all = blockSize(simulations@values, 
                     minblocks = 1, chunksize = chunksize)
  
  ## simulations@values in memory: bs_all$n = 1
  if (bs_all$n == 1){
    # extract data at locations
    y = matrix(simulations@values[locations,], byrow = TRUE, ncol = nPlumes(simulations))
    if (!regFun){
      interpolated = fun_interpol(y = y, data = newData[locations,], newdata = newData)  
    } else{
      interpolated = fun_interpol(y = y, data = newData[locations,], newdata = newData, dataLoc = locations, newdataLoc = 1:length(newData)) 
    }
    
    Interpolated = raster(interpolated,
                          xmn =  simulations@values@extent[1],
                          xmx =  simulations@values@extent[2],
                          ymn =  simulations@values@extent[3],
                          ymx =  simulations@values@extent[4],
                          crs =  crs(simulations@values))
  } else { # bs_all$n > 1; interpolated values do not fit into memory
    if (tmpfile == FALSE){
      stop("To process the data it needs to be saved to disk; 'tmpfile' must be a filename or path.")
    }
    #    if (bs_subset$n == 1){ # y does fit into memory
    #      Y = simulations@values[locations,, drop = FALSE] # does not work, produces wrong result
    #    } else { # bs_subset$n > 1; y does not fit into memory
    # save Y as raster
    Y = raster(nrow = length(locations),
               ncol = ncol(simulations@values),
               ext =  extent(simulations@values),
               crs =  crs(simulations@values))
    Y = writeStart(Y, 
                   filename = paste0(tmpfile, "_interpolateSubset.grd"),
                   overwrite = overwrite)
    h = 1
    for (i in 1:bs_all$n){
      # rows of this block 
      locations_i = bs_all$row[i] - 1 + 1:bs_all$nrows[i]
      which_locations_i = is.element(locations, locations_i)
      locations_i = locations[which_locations_i]
      
      locations_i = locations_i - (bs_all$row[i] - 1) 
      nLoc_i = length(locations_i)
      if (nLoc_i > 0){
        data_i = getValues(simulations@values, row = bs_all$row[i], nrows = bs_all$nrows[i]) 
        # indices to keep
        indices_i =  ncol(Y) * rep(locations_i - 1, each = ncol(Y)) + rep(1:ncol(Y), times = nLoc_i)
        subset_i = data_i[indices_i,]
        writeValues(Y, subset_i, start = h)
      }
      h = h + nLoc_i        
    }
    Y = writeStop(Y) 
    #}  
    # transpose to easier read data needed together
    y = t(Y)
    
    # raster to write interpolation result (transposed!)
    InterpolatedT = raster(nrow = ncol(simulations@values),
                           ncol = nrow(simulations@values),
                           ext =  extent(simulations@values),
                           crs =  crs(simulations@values))
    InterpolatedT = writeStart(InterpolatedT, 
                               filename = paste0(tmpfile, "_interpolateT.grd"),
                               overwrite = overwrite)
    
    bs_plumes = blockSize(InterpolatedT, 
                          minblocks = 1, chunksize = chunksize)
    
    for (i in 1:bs_plumes$n){
      y_i = matrix(getValues(y, row = bs_plumes$row[i], nrows = bs_plumes$nrows[i]),
                   byrow = FALSE, ncol = bs_plumes$nrows[i]) # transpose
    
      if (!regFun){
        interpolated_i = fun_interpol(y = y_i, data = newData[locations,], newdata = newData) 
      } else{
        interpolated_i = fun_interpol(y = y_i, data = newData[locations,], newdata = newData,
                                      dataLoc = locations, newdataLoc = 1:length(newData)) 
      }
      interpolated_i_numeric = as.numeric(interpolated_i)
      writeValues(InterpolatedT, interpolated_i_numeric, start = bs_plumes$row[i])
    }  
    InterpolatedT = writeStop(InterpolatedT)
    Interpolated = writeRaster(t(InterpolatedT), 
                               filename = paste0(tmpfile, "_interpolated.grd"),
                               overwrite = overwrite)
    rm(InterpolatedT)
    file.remove(paste0(tmpfile, "_interpolateT.grd"))
    file.remove(paste0(tmpfile, "_interpolateT.gri"))
    #if (bs_subset$n > 1){
    file.remove(paste0(tmpfile, "_interpolateSubset.grd"))
    file.remove(paste0(tmpfile, "_interpolateSubset.gri"))
    #}
  }
  return(Interpolated)
}


idw0z = replaceDefault(idw0, newDefaults = list(formula = z ~ 1))[["fun"]]
