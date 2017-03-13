length.SpatialDataFrame = function(x){nrow(x@data)} # should be superflous by proper definition
# idw0z = replaceDefault(idw0, newDefaults = list(formula = z ~ 1))[["fun"]] # defined in interpolate.R

interpol0reg = function(fun_interpolation = idw0z, # has to contain all parameters like formula, model... - may be a list
                               y, # matrix, must have same number of rows as data   matrix(simulations@values[locations,])                               
                               data, # SpatialDataFrame                             simulations@locations[locations,]
                               newdata, # SpatialDataFrame                          simulations@locations
                               dataLoc, # locations (indices) of given 'data' (relative to original 'simulations' it is taken from)
                               newdataLoc = 1:length(newdata),# locations (indices) of given 'newData' and 'y' (relative to original 'simulations' it is taken from)
                               dataSplit = 1:length(newdata), # list: indices or original 'simulations' to be used as input in this part
                               newdataSplit = 1:length(newdata)# list: indices of original 'simulations' to compute results for in this part
){
  parametersRaw = list()
  parList = rep(NA,3)

  
  # check class of input parameters
  # fun_interpolation
  switch(class(fun_interpolation),
         "function" = {
           parametersRaw[["fun_interpolation"]] = list()
           parametersRaw[["fun_interpolation"]][[1]] = fun_interpolation
         },
         "list" = {
           parametersRaw[["fun_interpolation"]] = fun_interpolation
           parList[1] = length(fun_interpolation)
         },
         {stop("'fun_interpolation' must be an (interpolation) function or a list of such functions.")}
  )
  # y
  if (!is.matrix(y)){
    stop("'y' must be a matrix (of input parameters).")
  }
  # data
  if (!is.SpatialDataFrame(data)){
    stop("'data' must be a 'SpatialDataFrame' with the locations of the input.")
  } else {
    if (nrow(y) != length(data) | nrow(y) != length(dataLoc)){
      stop("The rows in 'y', the locations in 'data' and 'dataLoc' refer to the same objects, 
            their numbers (", nrow(y),", ", length(data), ", ", length(dataLoc),") do not agree.")
    }
  }
  # newdata
  if (!is.SpatialDataFrame(newdata)){
    stop("'newdata' must be a 'SpatialDataFrame' with the locations for the output.")
  } else {
    if (length(newdata) != length(newdataLoc)){
      stop("The locations in 'newdata' and 'newdataLoc' refer to the same objects, 
            their numbers (", length(newdata), ", ", length(newdataLoc),") do not agree.")
    }
  }  
  # dataLoc
  if (!is.numeric(dataLoc)){
    stop("'dataLoc' has to be an integer vector of location indices.")
  }
  # newdataLoc
  if (!is.numeric(newdataLoc)){
    stop("'newdataLoc' has to be an integer vector of location indices.")
  }
  # dataSplit
  switch(class(dataSplit),
         "integer" = {
           parametersRaw[["dataSplit"]] = list()
           parametersRaw[["dataSplit"]][[1]] = dataSplit
         },
         "list" = {
           parametersRaw[["dataSplit"]] = dataSplit
           parList[2] = length(dataSplit)
         },
         {stop("'dataSplit' have to be integers (indices of locations) or a list of such.")}
  )
  # newdataSplit
  switch(class(newdataSplit),
         "integer" = {
           parametersRaw[["newdataSplit"]] = list()
           parametersRaw[["newdataSplit"]][[1]] = newdataSplit
         },
         "list" = {
           parametersRaw[["newdataSplit"]] = newdataSplit
           parList[3] = length(newdataSplit)
           if (any(duplicated(unlist(lapply(FUN = unique, X = newdataSplit))))){
             warning("The indices given in 'newdataSplit' overlap, in such locations the returned value belongs to the last given computation method.")
           }
         },
         {stop("'newdataSplit' have to be integers (indices of locations) or a list of such.")}
  )

  
  # compare length
  # are there any lists?
  if (all(is.na(parList))){ # no lists; turn all parameters into lists of length 1
    I = 1
  } else {
    # do all lists have same length?  
    I = unique(na.omit(parList))
    if (length(I) != 1){
      stop("The lists in 'model', 'beta', 'dataSplit', 'newdataSplit' do not all have the same length.")
    }
  }  
  
  
  # testing and preparation
  # are functions in fun_interpolation of correct type?
  for (i in seq(along = parametersRaw[["fun_interpolation"]])){
    funTest_i = replaceDefault(fun = parametersRaw[["fun_interpolation"]][[i]], type = "fun_interpolation.interpolate")
    if (!funTest_i[["accept"]]){
      stop(paste0("The ", i, "-th of the functions in 'fun_interpolation' is not an interpolation function."))
    } else {
      parametersRaw[["fun_interpolation"]][[i]] = funTest_i[["fun"]]
    }
  }

  # which rows of data / newdata should be used in which interpolation?
  # those that are part of the respective dataSplit / newdataSplit
  for (i in seq(along = parametersRaw[["dataSplit"]])){
    parametersRaw[["dataSplit"]][[i]] = which(is.element(dataLoc, parametersRaw[["dataSplit"]][[i]]))
  }
  for (i in seq(along = parametersRaw[["newdataSplit"]])){
    parametersRaw[["newdataSplit"]][[i]] = which(is.element(newdataLoc, parametersRaw[["newdataSplit"]][[i]]))
  }
  
  # if single elements were given instead of lists, make copies
  parameters = list()
  for (j in 1:3){
    parameters[[j]] = list()
    if (is.na(parList[j])){ # parameter is not a list -> make list of required number of copies
      for (k in 1:I){
        parameters[[j]][[k]] = parametersRaw[[j]][[1]]  
      }
    } else { # parameter is a list -> take it as it is
      parameters[[j]]  = parametersRaw[[j]] 
    }
  }
  names(parameters) = names(parametersRaw)

  # result must be matrix, nrow like newdata, ncol like y
  result = matrix(nrow = length(newdata), ncol = ncol(y))
  
  # interpolate for all list entries
  for (i in 1:I){
    # result[parameters[["newdataSplit"]], ] = 
    #          fun_interpol(formula = formula,
    #                       data = data[parameters[["dataSplit"]][[i]],],
    #                       newdata = newdata[parameters[["newdataSplit"]][[i]],],
    #                       y = y[parameters[["dataSplit"]][[i]],],
    #                       model = parameters[["model"]][[i]])
    #debug(parameters[["fun_interpolation"]][[i]])
    result[parameters[["newdataSplit"]][[i]], ] = 
      parameters[["fun_interpolation"]][[i]](data = subsetSDF(data, locations = parameters[["dataSplit"]][[i]]),
                           newdata = subsetSDF(newdata, locations = parameters[["newdataSplit"]][[i]]),
                           y = y[parameters[["dataSplit"]][[i]],,drop = FALSE])
  }
  return(result)
}
