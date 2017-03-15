############################################################
#        extractSpatialDataFrame                           #
############################################################
# extractSpatialDataFrame = function(
#   obj,
#   kinds = 1:nlayers(obj@values), # integer, layers to be returned 
#   plumes = 1:ncol(obj@values), # integer, plumes to be returned
#   chunksize = 1e+7
# ){}

extractSpatialDataFrame = function(
  obj, # Simulations object
  kinds = 1:nlayers(obj@values), # integer, layers to be returned 
  plumes = 1:ncol(obj@values), # integer, plumes to be returned
  chunksize = 1e+7
){
  #library(raster)
  nL = nrow(obj@values)
  
  nK = nlayers(obj@values)
  kindLength0 = length(kinds)
  kinds = intersect(kinds, 1:nK)
  kindLength1 = length(kinds)
  if (kindLength1 < kindLength0){
    warning("'kinds' contains invalid numbers, they are ignored.")
    if (kindLength1 == 0){
      stop("No valid 'kinds'.")
    }
  }
  
  nP = ncol(obj@values)
  plumeLength0 = length(plumes)
  plumes = intersect(plumes, 1:nP)
  plumeLength1 = length(plumes)
  if (plumeLength1 < plumeLength0){
    warning("'plumes' contains invalid numbers, they are ignored.")
    if (plumeLength1 == 0){
      stop("No valid 'plumes'.")
    }
  }
  if (nL > chunksize){
    stop("Already one map (one kind of one plume) is too big to fit into memory.")
  }
  ncellsAll = kindLength1 * plumeLength1 * nL
  if (ncellsAll > chunksize){ 
    newKindLength = max(floor(kindLength1/ncellsAll * chunksize), 1)
    ncellsAll = newKindLength * plumeLength1 * nL
    if (ncellsAll > chunksize){
      newPlumeLength = max(floor(plumeLength1/ncellsAll * chunksize), 1)
      ncellsAll = newKindLength * newPlumeLength * nL
      if (ncellsAll > chunksize){
        stop(paste0("chunksize: ", chunksize, "; newKindLength: ", newKindLength, "; newPlumeLength: ", newPlumeLength))
      }  
      if (newPlumeLength < plumeLength1){
        warning(paste0("Only the first ", newPlumeLength, " indices of 'plumes' are used."))
        plumes = plumes[1:newPlumeLength]  
      }      
    }
    if (newKindLength < kindLength1){
      warning(paste0("Only the first ", newKindLength, " indices of 'kinds' are used."))
      kinds = kinds[1:newKindLength]
    }    
  }
  
  kindLength = length(kinds)
  plumeLength = length(plumes)
  valuesLayer = subset(obj@values, kinds)
  valuesAllPlumes = getValues(valuesLayer)
  if (kindLength == 1){
    valuesAllPlumes = matrix(valuesAllPlumes, ncol = 1)
  }
  values = valuesAllPlumes[
    rep(plumes, times = nL) + nPlumes(obj) * rep(0:(nL-1), each = plumeLength),, drop = FALSE]
  
  Values = matrix(nrow = nL, ncol = plumeLength * kindLength)
  for (i in 1:kindLength){
    Values[,(i - 1) * plumeLength + 1:plumeLength] = 
      matrix(values[,i], nrow = nL, byrow = TRUE)
  }
  objLoc = obj@locations
  objLoc@data = as.data.frame(Values)
  names(objLoc@data) = paste0("plume", rep(plumes, times = kindLength), 
                              "_kind", rep(kinds, each = plumeLength))
  return(objLoc)
}  
 
