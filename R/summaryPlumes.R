#########################################################
#                    summaryPlumes                      #
#########################################################


weightedMean = function(x, weight, na.rm = FALSE
){
  weightedMean = mean(x * weight, na.rm = na.rm)
  return(weightedMean)
}


# decision: no user-defined functions: associativeness is a strong restriction, also the existence of a neutral element
## for really free user-defined functions (e.g. multiple detection at a certain threshold) use simulationsApply
summaryPlumes = function(
  simulations,
  locations = 1:nLocations(simulations),
  plumes = 1:nPlumes(simulations),
  kinds = 1,
  fun, # must be a function that can be applied to vectors and is associative fun(c(a,b)) = fun(c(fun(a),b)), it has to return one single value
  summaryFun = weightedMean,
  weight = 1, # either a vector of same length as plumes (i.e. for the plumes to actually be used) or a character indicating a column of simulations@plumes
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

  #  select values of interest
  if (length(kinds) > 1){
    warning("Only the first of the 'kinds'(", names(data)[kinds[1]], ") was used.", sep = "")
    kinds = kinds[1]
  }
  data = data[[kinds]]
#   if (!missing(kinds)){
#     data = data[[kinds[1]]]
#     if (length(kinds) > 1){
#       warning("Only the first of the 'kinds'(", names(data)[kinds[1]], ") was used.", sep = "")
#     }
#   } else {
#     warning("No 'kinds' was indicated, therefore the first kind (", names(data)[1], ") was used.", sep = "")
#     data = data[[1]]
#   }

  # determine plumes to be used
  nP = ncol(data)
  if (!identical(plumes, 1:nP)){
    plumesIn = plumes > 0 & plumes <= nP
    plumes[!plumesIn] = NA
    if (any(is.na(plumes))){
      warning("Some of the selected 'plumes' are out of bounds or 'NA', these values are set to 'NA'.")
    }
  }
  nPl = length(plumes)

  if (class(simulations) == "Simulations"){
    if (class(weight) == "character"){
      weight = simulations@plumes[plumes, weight]
    }
  }

  # determine locations to be used
  nL = nrow(data)
  isLocations = FALSE
  nLocNA = 0

  if (!identical(locations, 1:nL)){
    locationsIn = locations > 0 & locations <= nL
    locations[!locationsIn] = NA
    if (any(is.na(locations))){
      warning("Some of the selected 'locations' are out of bounds or 'NA', these values are set to 'NA'.")
    }
    nLocNA = sum(is.na(locations))
    locations = sort(locations)
    nL = length(locations)
    isLocations = TRUE
  }

  # block-wise processing
  data0 = matrix(nrow = 0, ncol = nPl)
  summary = apply(FUN = fun, X = data0, MARGIN = 2)

  bs = blockSize(data, minblocks = 1, n = 1)
  message(paste("Data is processed in", bs$n, "block(s)."))

  # select the data to actually be used
  for (j in 1:bs$n) {
    theseLoc = bs$row[j] - 1 + 1:bs$nrows[j]
    if (isLocations){
      theseLoc = intersect(locations, theseLoc) #- (bs$row[j] - 1)
    }
    nLoc = length(theseLoc)

    if(nLoc >= 1){ # if no rows from this block are used, jump to next block
      in_j = getValues(data, row = bs$row[j], nrows = bs$nrows[j])
      index = selectRow(theseLoc - bs$row[j] + 1, bs$nrows[j], nP)
      index = index[selectCol(plumes, nLoc, nP)]

      selected_j = in_j[index]
      summary = apply(FUN = fun, MARGIN = 2, na.rm = na.rm,
                      X = matrix(c(summary, selected_j), ncol = nPl, byrow = TRUE))
    }
  }
  if (nLocNA > 0){
    summary = apply(FUN = fun, MARGIN = 2, na.rm = na.rm,
                    X = matrix(c(summary, rep(NA, nPl * nLocNA)), ncol = nPl, byrow = TRUE))
  }
  # global summary, using weights
  summary[is.na(plumes)] = NA
  globalSummary = summaryFun(summary, weight = weight)
  out = list()
  out[["summary"]] = globalSummary
  names(summary) = plumes
  out[["summaryPlumes"]] = summary
  return(out)
}
