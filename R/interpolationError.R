###############################################################################
# interpolationError                                                               #
###############################################################################


# interpolation error: interpolate, apply error function via simulationsApply
interpolationError = function(
  simulations,
  locations,
  kinds,
  fun_interpolation = NA,
  fun_error = NA,
  fun_Rpl = NA,
  fun_Rpl_cellStats = "mean",
  fun_l = NA,
  tmpfile = "tmp_interpolationError",
  overwrite = FALSE,
  chunksize = 1e+7
  ){
  # keep only relevant layer of simulations@values
  simulationsSubset = simulations
  simulationsSubset@values = subset(simulations@values, kinds[1], drop = FALSE)
  
  # check functions (except those that are tested inside further functions)
  useRpl = FALSE
  if (!is.na(fun_Rpl)){
    RplTested = replaceDefault(fun_Rpl, type = "fun.simulationsApply")
    bs_fun_pl = blockSize(simulationsSubset@values, 
                          n = formals(fun_error)[["nout"]], 
                          minblocks = 1, chunksize = chunksize)
    if (!RplTested[["accept"]] | bs_fun_pl$n > 1) {
      warning("'fun_Rpl' invalid or cannot be applied as the result of 'fun_error' does not fit into memory, 
              the value returned as 'cost' is the result of 'fun_Rpl_cellStats'.")
    } else {
      useRpl = TRUE
    }
  }
  useL = FALSE
  if (!is.na(fun_l)){
    lTested = replaceDefault(fun_l, type = "fun.simulationsApply")
    if (!lTested[["accept"]]) {
      warning("'fun_l' invalid, no result returned.")
    } else {
      useL = TRUE
    }
  }  
  # check class of simulations
  if (!is(simulations, "Simulations")){
    stop("'simulations' must be of class 'Simulations'.")
  }
 
  # interpolate
  interpolated = interpolate(
    simulations = simulationsSubset,
    locations = locations,
    fun_interpolation = fun_interpolation,
    tmpfile = tmpfile,
    overwrite = overwrite,
    chunksize = chunksize
  )
    
  # generate simulations of original and interpolated values
  simulationsSubset@values = stack(simulationsSubset@values,
                                   interpolated)
  
  # apply error functions value-wise and summary function(s) via simulationsApply
  if (useRpl){
    if (useL){
      error = simulationsApply(
        simulations = simulationsSubset,
        fun_pl = fun_error,
        fun_Rpl = fun_Rpl,
        fun_Rpl_cellStats = fun_Rpl_cellStats,
        fun_l = fun_l,
        tmpfile = tmpfile, 
        overwrite = overwrite,
        chunksize = chunksize
      )       
    } else {
      error = simulationsApply(
        simulations = simulationsSubset,
        fun_pl = fun_error,
        fun_Rpl = fun_Rpl,
        fun_Rpl_cellStats = fun_Rpl_cellStats,
        tmpfile = tmpfile, 
        overwrite = overwrite,
        chunksize = chunksize
      )  
    }  
  } else {
    if (useL){
      error = simulationsApply(
        simulations = simulationsSubset,
        fun_pl = fun_error,
        fun_Rpl_cellStats = fun_Rpl_cellStats,
        fun_l = fun_l,
        tmpfile = tmpfile, 
        overwrite = overwrite,
        chunksize = chunksize
      )      
    } else {
      error = simulationsApply(
        simulations = simulationsSubset,
        fun_pl = fun_error,
        fun_Rpl_cellStats = fun_Rpl_cellStats,
        tmpfile = tmpfile, 
        overwrite = overwrite,
        chunksize = chunksize
      )     
    }
  }

  
  # put together result
  result = list()
  if (useRpl){
    result[["cost"]] = error[["result_global_locationsplumes"]]
    result[["cost_cellStats"]] = error[["cellStats_global_locationsplumes"]]
  } else {
    result[["cost"]] = error[["cellStats_global_locationsplumes"]]
  }
  result[["error_locationsplumes"]] = error[["result_locationsplumes"]]
  result[["interpolated"]] = interpolated
  if (useL){
    result[["costLocations"]] = error[["result_locations"]]
  }

  return(result) 
}
  
