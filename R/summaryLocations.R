#################################################################
#                summaryLocations                               #
#################################################################

summaryLocations = function(
  simulations,
  locations = 1:nLocations(simulations), # does not work for simulations of class raster
  plumes = 1:nPlumes(simulations), # does not work for simulations of class raster
  kinds = 1,
  fun, # must be a function that can be applied to vectors, it has to return one single value
  summaryFun = weightedMean,
  weight = 1, # either a vector of same length as locations (i.e. for the locations to actually be used) or a character indicating a column of simulations@locations
  na.rm = FALSE,
  ...
){
  # adapt to class of simulations
  if (is.element(class(simulations), c("RasterLayer", "RasterStack", "RasterBrick"))){
    data = simulations
  } else {
    if (class(simulations) == "Simulations"){
      data = simulations@values
    } else {
      stop("'simulations' must be of class 'Simulations' or of a 'Raster*' class.")
    }
  }
  
  # select values of interest
  #if (!missing(kinds)){
  data = subset(data, kinds[1])
  if (length(kinds) > 1){
    warning("Only the first of the 'kinds'(", names(data)[kinds[1]], ") was used.", sep = "")
  }
  #} 
  #! below not necessary, just uses default
  #else {
  #  warning("No 'values' was indicated, therefore the first kind (", names(data)[1], ") was used.", sep = "")
  #  data = subset(data, values[1])
  #}
  layerNames = names(data)
  
    
  # determine plumes to be used
  nP = ncol(data)  
  #if (!missing(plumes)){
    plumes[is.na(plumes)] = 0
    plumesIn = plumes > 0 & plumes <= nP
    if (any(!plumesIn)){
      warning("Some of the selected 'plumes' are out of bounds or 'NA', these values are set to 'NA'.")
    }
    plumes[!plumesIn] = NA
    nPl = length(plumes)
  #}else{
  #  plumes = 1:nP
  #  nPl = nP
  #}

  # determine locations to be used
  nL = nrow(data)
  #isLocations = FALSE
  nLocNA = 0
  #if (!missing(locations)){
    locations[is.na(locations)] = 0
    locationsIn = locations > 0 & locations <= nL
    if (any(!locationsIn)){
      warning("Some of the selected 'locations' are out of bounds or 'NA', these values are set to 'NA'.")
    }
    locations[!locationsIn] = NA
    if (class(simulations) == "Simulations"){
      if (is.character(weight)){
        weight = simulations@locations@data[locations,weight]
      }
    } 
    locationsTable = table(locations)
    locationsRank = rank(locations)
    nLocNA = sum(!locationsIn)
    locations = sort(unique(locations))
    nL = length(locations)
    if (identical(locations, 1:nL)){
      isLocations = FALSE
    } else {
      isLocations = TRUE 
    }
  #} 
  if (!isLocations & class(simulations) == "Simulations" & is.character(weight)){
    weight = simulations@locations@data[,weight]
  }
  
  # block-wise processing
  summary = rep(NA, nL)
  
  bs = blockSize(data, minblocks = 1, n = 1)
  message(paste("Data is processed in", bs$n, "blocks."))
  
  # select the data to actually be used 
  k = 0
  for (j in 1:bs$n) {
    theseLoc = bs$row[j] - 1 + 1:bs$nrows[j]
    if (isLocations){
      theseLoc = intersect(locations, theseLoc)# - (bs$row[j] - 1) 
    }
    nLoc = length(theseLoc)    
    
    if(nLoc >= 1){ # if no rows from this block are used, jump to next block
      in_j = getValues(data, row = bs$row[j], nrows = bs$nrows[j])
      index = selectRow(theseLoc - bs$row[j] + 1, bs$nrows[j], nP)
      index = index[selectCol(plumes, nLoc, nP)]
      
      selected_j = in_j[index]
      summary[k + 1:nLoc] = apply(FUN = fun, MARGIN = 1, na.rm = na.rm,
                      X = matrix(selected_j, ncol = nPl, byrow = TRUE))
      k = k + nLoc
    }     
  }
  # restore original order and multiplicity
  if(isLocations){
    names(summary) = locations
    summary = unlist(mapply(rep, summary, locationsTable))[locationsRank]    
  }

  # global summary, using weights  
  globalSummary = summaryFun(summary, weight = weight)
  out = list()
  out[["summary"]] = globalSummary  
  out[["summaryLocations"]] = summary
  return(out)
}
