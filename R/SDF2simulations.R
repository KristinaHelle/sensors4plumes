#SDF2simulations = function(x, 
#                           indices = matrix(1:ncol(x@data), nrow = 1)){}

SDF2simulations = function(
  x, # SpatialDataFrame
  indices = matrix(1:ncol(x@data), nrow = 1) # indices of columns of x to be transformed into values, 
  #other values are kept in locations; can be matrix, then each row contains the data for one layer
  ){
  if (!is.matrix(indices)){
    stop("'indices' must be a matrix, each row indicating the columns that belong to the same kind of values.")
  }
  if (any(!is.element(as.integer(indices), 1:ncol(x@data)))){
    stop(paste0("'indices' must only contain values between 1 and ", ncol(x@data), " and no 'NA'."))
  }
  noIndices = setdiff(1:ncol(x@data), indices)
  Values = list()
  for (i in 1:nrow(indices)){
    Values[[i]] = raster(as.matrix(x@data[indices[i,]]),
                         xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                         crs = "+init=epsg:4326")
  }
  names(Values) = row.names(indices)
  Plumes = as.data.frame(t(indices))
  Locations = x
  Locations@data = Locations@data[, noIndices]
  out = Simulations(
    locations = Locations,
    plumes = Plumes,
    values = stack(Values)
  )
  # data type of all marked columns must fit (else adapt to most complex)
  # if toValues is matrix, rows need to have same number of entries
  # names of columns are kept in "plumes"
  # layer names can be taken from rownames of 'toValues'
  return(out)
}
#setMethod("SDF2simulations", signature(x = "SpatialDataFrame"), SDF2simulations.SpatialDataFrame)