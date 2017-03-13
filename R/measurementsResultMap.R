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

singleDetection = replaceDefault(measurementsResult, newDefaults = list(
  kinds = "detectable",
  fun_p = function(x, nout = 1)
    {
    prod(1-x)
    },
  fun_Rp = function(x, weight = 1, nout = 1)
    {
    mean(x * weight$totalDose)/mean(weight$totalDose)
  },
  fun_pl = function(x, nout = 1){
    x
  },
  fun_Rpl = function(x, weight_l = 1, weight_p = 1, nout = nLocations(simulations))
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
  }),
  type = "costFun.optimiseSD")[[1]]  

multipleDetection = replaceDefault(measurementsResult, newDefaults = list(
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
  type = "costFun.optimiseSD")[[1]] 

earlyDetection = replaceDefault(measurementsResult, newDefaults = list(
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
  type = "costFun.optimiseSD")[[1]]

