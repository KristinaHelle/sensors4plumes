data(radioactivePlumes)

# compute and add often used properties
## values
threshold = 1e-7
radioactivePlumes@values$detectable = calc(
  radioactivePlumes@values$maxdose,
  fun = function(x){x >= threshold})

## locations
radioactivePlumes@locations@data$area = as.numeric(table(radioactivePlumes@locations@index)/16)
radioactivePlumes@locations@data$index = NULL

## plumes
names(radioactivePlumes@plumes)[1] = "date"
sumWeightArea = function(x, weight = radioactivePlumes@locations@data$area, nout = 1){
  sum(x * weight)
}
radioactivePlumes@plumes$totalDose = simulationsApply(simulations = radioactivePlumes,
                                                      fun_p = sumWeightArea, kinds = "finaldose")[["result_plumes"]][,1]
radioactivePlumes@plumes$nDetectable = 
  summaryPlumes(radioactivePlumes, fun = sum, kinds = "detectable")[["summaryPlumes"]]

