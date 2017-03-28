##############################################################
#  measurementResult                                         #
##############################################################
# cost functions based on 'measurements' only and via plume-wise cost (detection)

# basic function: measurement result
# specifications:
## singleDetection
## multipleDetection
## earlyDetection
measurementsResult = function(
  simulations,
  locations,
  kinds,
  fun_p = NA,
  fun_Rp = NA,
  fun_pl = NA,
  fun_Rpl = NA
#   fun_p = NA,
#   fun_Rp = NA,
#   fun_pl = NA,
#   fun_Rpl = NA
){
  result = simulationsApply(
    simulations = simulations,
    locations = locations,
    kinds = kinds,
    fun_p = fun_p,
    fun_Rp = fun_Rp
  )
  out = list()
  out[["cost"]] = result[["result_global_plumes"]]
  if (is.function(fun_pl)){
    simulations@plumes$result_plumes =
      result[["result_plumes"]]

    map = simulationsApply(
      simulations = simulations,
      kinds = kinds,
      fun_pl = fun_pl,
      fun_Rpl = fun_Rpl
      )

  out[["costLocations"]] = map[["result_global_locationsplumes"]]
  }
  out[["costPlumes"]] = result[["result_plumes"]]
  return(out)
}


singleDetection = function(simulations, locations, plot = FALSE){
 x1 = function(x, nout = 1){
   x
 }
 prodNeg1 = function(x, nout = 1){
   prod(1-x)
 }
 meanWeight_totalDose1 = function(x, weight = 1, nout = 1)
 {
   mean(x * weight$totalDose)/mean(weight$totalDose)
 }
 sumUndetected = function(x, weight_l = 1, weight_p = 1, nout = nLocations(simulations))
 {
   nout = nrow(weight_l)
   detectableAsMatrix = matrix(x, nrow = nrow(weight_l), byrow = TRUE)
   detectableUndetected = detectableAsMatrix[,which(weight_p$result_plumes == 1), drop = FALSE]
   if (dim(detectableUndetected)[2] >= 1){
     map = apply(FUN = sum, X = detectableUndetected, MARGIN = 1)
   } else {
     map = rep(0, nout)
   }

   return(map)
 }
 if (plot){
   sD = replaceDefault(measurementsResult, newDefaults = list(
     kinds = "detectable",
     fun_p = prodNeg1,
     fun_Rp = meanWeight_totalDose1,
     fun_pl = x1,
     fun_Rpl = sumUndetected),
     type = "costFun.optimiseSD")[[1]](simulations = simulations, locations = locations)
 } else {
   sD = replaceDefault(measurementsResult, newDefaults = list(
     kinds = "detectable",
     fun_p = prodNeg1,
     fun_Rp = meanWeight_totalDose1),
     type = "costFun.optimiseSD")[[1]](simulations = simulations, locations = locations)
 }
  return(sD)
}

multipleDetection = function(simulations, locations, plot = FALSE){
  if (plot){
    warning("'plot = FALSE' is set as no plotting function implemented.")
  }
    mD = replaceDefault(measurementsResult, newDefaults = list(
      kinds = "detectable",
      fun_p = function(x, nout = 1)
      {
        sum(x)
      },
      fun_Rp = function(x, weight = 1, nout = 1)
      {
        y = (weight$nDetectable - x)/weight$nDetectable
        mean(y)
      }),
      type = "costFun.optimiseSD")[[1]](simulations = simulations, locations = locations)

  return(mD)
}

earlyDetection = function(simulations, locations, plot = FALSE){
  if (plot){
    warning("'plot = FALSE' is set as no plotting function implemented.")
  }
  eD = replaceDefault(measurementsResult, newDefaults = list(
    kinds = "time",
    fun_p = function(x, nout = 1)
    {
      minX = min(x, na.rm = TRUE)
    },
    fun_Rp = function(x, nout = 1, weight = 1,
                      notDetected = 6.0480e+05 * 2)
    {
      xFinite = x - weight$earliestDetection
      xFinite[is.infinite(x)] = notDetected
      out = mean(xFinite)/notDetected
    }),
    type = "costFun.optimiseSD")[[1]](simulations = simulations, locations = locations)
  return(eD)
}
