#############################################################
#           replaceDefault                                  #
#############################################################

# parameters
# - Parameters: parameters of the function
# - Required: needed (with or without default (default will be overwritten)), depends on type
# - RequiredDefault: needed with default ('nout' in fun.simulationsApply - to use outside the function), depends on type
# - DefaultOld: given as old default
# - DefaultNew: given as new default
# 
# by definition
# DefaultOld SUBSET Parameters
# RequiredDefault SUBSET Required
# 
# test
# (Required AND RequiredDefault) SUBSET Parameters   # all required parameters exist?
# DefaultNew SUBSET Parameters                       # all newDefaults can be used as parameters?
# RequiredDefault SUBSET (DefaultOld AND DefaultNew) # all required defaults set?


replaceDefault = function(fun, newDefaults = list(), type){
  #warn = character(0)
  acceptFun = TRUE
  if (!is.function(fun)) {   
    #warn = c(warn, "must be a function.")
    warning("must be a function")
    acceptFun = FALSE
  } else {
    # names of existing/required parameters and defaults
    parameters = names(formals(fun))
    oldDefault = parameters[!sapply(formals(fun), is.name)] # parameters with default values
    newDefault = names(newDefaults)
    if (!missing(type)){
      required = switch(type,
                        "summaryFun.summaryPlumes" = c("x", "weight"),
                        "fun.simulationsApply" = c("x", "nout"),
                        "funR.simulationsApply" = c("x", "nout", "weight"),
                        "funRR.simulationsApply" = c("x", "nout", "weight_l", "weight_p"),
                        "fun_interpolation.interpolate" = c("data", "newdata", "y"),
                        #"fun_interpolationSplit.interpolate" = c("data", "newdata", "y", "dataLoc", "newdataLoc"),
                        "fun.spatialSpread" = c("allLocations", "locations"),
                        "fun_R.spatialSpread" = c("x"),
                        "costFun.optimiseSD" = c("simulations", "locations"),
                       # "costMap.optimiseSD" = c("simulations", "locations", "nameSave = nameSave"), # in use?
                        "optimisationFun.optimiseSD" = 
                          c("simulations", "costFun", 
                            "locationsAll", "locationsFix", "locationsInitial",
                            "aimCost", "aimNumber", "nameSave"
                          ),
                        "evalFunc.rbga.bin" = c("chromosome"),
                        #"rbga.bin" = c("size", "suggestions", "popSize", "iters", "mutationChance", "elitism", "zeroToOneRatio", "monitorFunc", "evalFunc", "showSettings", "verbose"),
                        character(0)
      )        
      requiredDefault = setdiff(parameters, required) # all not-required parameters need default values
      
      if (type == "fun.simulationsApply" | type == "funR.simulationsApply"){
        requiredDefault = c(requiredDefault, "nout")
      }   
    } else {
      required = character(0)
      requiredDefault = character(0)
    }
    if (!missing(type)){
      # do parameters cover all required values?
      requiredAreParameters = is.element(required, parameters)
      if (!all(requiredAreParameters)){
        #warn = c(warn, paste0(c("function parameters missing:", required[!requiredAreParameters]), collapse = " "))
        warning(paste0(c("function parameter(s) missing:", required[!requiredAreParameters]), collapse = " "))
        acceptFun = FALSE
      }
      # are all required defaults given?
      requiredDefaultsDefined = is.element(requiredDefault, union(oldDefault, newDefault))
      if (!all(requiredDefaultsDefined)){
        #warn = c(warn, paste0(c("default values missing for the parameters:", requiredDefault[!requiredDefaultsDefined])))
        warning(paste0(c("default values missing for the parameter(s):", requiredDefault[!requiredDefaultsDefined]), collapse = " "))
        acceptFun = FALSE
      }
      # are default values given for non-parameters?
      defaultHasParameter = is.element(newDefault, parameters)
      if (!all(defaultHasParameter)){
        #warn = c(warn, paste0(c("in 'newDefault' values given for non-existing parameters:", newDefaults[!defaultHasParameter])))
        warning(paste0(c("in 'newDefaults' values given for non-existing parameter(s):", 
                         newDefault[!defaultHasParameter]), collapse = " "))
        acceptFun = FALSE
      }
    }

    parameterHasNewDefault = is.element(parameters, newDefault)
    for (i in seq(along = parameters)){
      if (parameterHasNewDefault[i]){
        formals(fun)[[parameters[i]]] = newDefaults[[parameters[i]]]
      }
    }    
  }
  out = list(fun = fun, accept = acceptFun)#, warnings = warn)
#  if (length(warn) == 0){
#    return(fun)    
#  } else {
#    return(warn)
#  }
  return(out)
}     
