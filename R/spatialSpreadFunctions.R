########################################################
# spatial spread functions                             #
########################################################
# minimalDistance
# krigingVariance
# to be used as 'fun' in 'spatialSpread'
## it would be nice to have more like Entropy, but algorithm not known



minimalDistance = function(allLocations, 
                           locations, 
                           algorithm = "kd_tree"){
  nLoc = length(locations)
  coordinatesAll = coordinates(allLocations)
  nAll = nrow(coordinatesAll)
  
  if (nLoc == 0) {
    minDist = rep(Inf, nAll)
  }
  if (nLoc == 1) {
    minDist = sqrt((coordinatesAll[,1] - coordinatesAll[locations,1]) ^ 2 + 
                     (coordinatesAll[,2] - coordinatesAll[locations,2]) ^ 2)
  }  
  if (nLoc > 1) {
    minDist = get.knnx(data = coordinatesAll[locations,], 
                       query = coordinatesAll, k = 1, 
                       algorithm = algorithm)$nn.dist[,1]              
  }
  
  result = list()
  result[["cost"]] = minDist
  return(result)
}



krigingVariance = function(allLocations,
                           locations,
                           model # must be set by replaceDefault
) {
  nLoc = length(locations)
  coordinatesAll = coordinates(allLocations)
  nAll = nrow(coordinatesAll)
  
  if (nLoc == 0){
    krigVar = rep(Inf, nAll)
  } else {
    # make dummy data for kriging
    allLocations@data$z = 1
    krigVar = krige(z ~ 1, 
                    allLocations[locations,,drop = FALSE],
                    allLocations,
                    model = model)$var1.var
  }
  result = list()
  result[["cost"]] = krigVar
  return(result)
}
