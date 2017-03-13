############################################################
# locationsCost                                            #
############################################################
# cost that is based only on the coordinates of the locations
## wrapper function to clean input


spatialSpread = function(
  simulations,
  locations, 
  weightByArea = TRUE,
  fun = NA,
  fun_R = NA 
  ){
  # test function
  fun_ = replaceDefault(fun, type = "fun.spatialSpread")
  if (!fun_[["accept"]]){
    stop("'fun' does not have the correct parameters.")
  }
  useFunR = FALSE
  if (is.function(fun_R)){
    fun_R_ = replaceDefault(fun_R, type = "fun_R.spatialSpread")
    if (!fun_R_[["accept"]]){
      warning("'fun_R' does not have the correct parameters, it is not applied.")     
    } else {
      useFunR = TRUE
    }   
  }

  # get coordinates from simulations
  ## simulations can be "Simulations" or "SpatialDataFrame" 
  if (is(simulations, "Simulations")){
    all = simulations@locations
  } else {
    if (is(simulations, "SpatialDataFrame")){
      all = simulations
    } else {
      stop("'simulations' must be of class 'Simulations' or 'SpatialDataFrame'.")
    }
  }
  ## all must be "Spatial"
  if (is(all, "SpatialIndexDataFrame")){
    stop("The locations of 'simulations' are a 'SpatialIndexDataFrame', 
         no coordinates-based cost can be computed.")
  }
  areas =  areaSDF(all)
  if (is(all, "SpatialPolygridDataFrame")){
    all = as(all, "SpatialPointsDataFrame")
  }
  
  # clean locations
  locations[locations < 1 | locations > nrow(all@data)] = NA
  locations = unique(locations[!is.na(locations)])
  
  # compute location-wise cost
  cost = fun_[["fun"]](allLocations = all,
               locations = locations)
  
  # (summarise and) return 
  result = list()
  if (useFunR) {
    if (weightByArea){   
      locationwiseCost = cost[["cost"]] * areas
    } else {
      locationwiseCost = cost[["cost"]]
    }
    result[["cost"]] = fun_R_[["fun"]](x = locationwiseCost)    
  }
  result[["costLocations"]] = cost[["cost"]]
  return(result)
}  
 

