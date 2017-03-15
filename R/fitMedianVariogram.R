##########################################################
# fitMedianVariogram                                     #
##########################################################
# needs automap: autoFitVariogram
fitMedianVariogram = function(simulations, plumes, locations, kinds = 1){
  if (!missing(locations)){
    # clean locations
    locationsP = locations
    locations[locations < 1 | locations > nLocations(simulations)] = NA
    locations = locations[!is.na(locations)]
    locations = sort(unique(locations))
#    locUnique = FALSE
#    if (length(locations) == length(locationsP)){
#      locUnique = TRUE
#    }    
#  } else {
#    locUnique = TRUE
  }

  # subset data and turn into Spatial [sp] object
  if (!missing(locations)){
    samplePlumesValues0 = subset(simulations, kinds = kinds[1], 
                                 locations = locations, plumes = plumes)    
  } else {
    samplePlumesValues0 = subset(simulations, kinds = kinds[1], 
                                 plumes = plumes)
  }

  samplePlumesValues = extractSpatialDataFrame(samplePlumesValues0)
  if (is (samplePlumesValues, "SpatialPolygridDataFrame")){
#    if (locUnique){
#      samplePlumesValues = as(samplePlumesValues, "SpatialGridDataFrame")      
#    } else {
      samplePlumesValues = as(samplePlumesValues, "SpatialPointsDataFrame")
#    }

  }
  if (is (samplePlumesValues, "SpatialIndexDataFrame")){
    stop("Variogram cannot be computed as locations of 'simulations' are 'SpatialIndexDataFrame'.")
  }
#  samplePlumesValues = as(samplePlumesValues1, "SpatialPointsDataFrame")

  # fit variograms
  variograms = list()
  for (i in seq(along = plumes)) {
#    if (locUnique){
#      thisPlume = samplePlumesValues[,,i]      
#    } else { 
      thisPlume = samplePlumesValues[,i]
#    }

    thisFormula = as.formula(paste(names(thisPlume@data), "~ 1"))
    variograms[[i]] = autofitVariogram(formula = thisFormula, input_data = thisPlume)
  }
  
  # extract parameters
  variogramParameters  = array(dim = c(2, 9, length(plumes)))
  dimnames(variogramParameters)[[2]] = 
    c("model", "psill", "range", "kappa", "ang1", "ang2", "ang3", "anis1", "anis2")
  for (i in seq(along = plumes)){
    parametersHere = names(variograms[[i]][[2]])
    for (j in seq(along = parametersHere)){
      variogramParameters[,parametersHere[j],i] = 
        variograms[[i]][[2]][,parametersHere[j]]
    } 
  }
  # generate median variogram
  variogramParametersSummary = apply(FUN = summary, X = variogramParameters, MARGIN = c(1,2))
  
  medianVariogram = vgm(psill = variogramParametersSummary[3,2,"psill"],
                      model = levels(variograms[[1]]$var_model$model)[
                        modal(variogramParameters[2,"model",])],
                      range = variogramParametersSummary[3,2,"range"],
                      nugget = variogramParametersSummary[3,1,"psill"],
                      kappa = variogramParametersSummary[3,2,"kappa"])
  return(medianVariogram)
}





