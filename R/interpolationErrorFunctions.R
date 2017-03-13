###########################################################
# interpolation error functions                           #
###########################################################

# functions to compare original and map point-wise (to be applied via simulationsApply)

## pl: point-and-plume-wise
# absolute error
absError = function(x, nout = 1){
  out = abs(x[1] - x[2])
}

# delineation error (version that keeps 3 values can only be used as cost function if fun_Rpl only takes one layer)
delineationError = function(x, nout = 1,#nout = 3, 
                            threshold = 1e-7, weightFalseNeg = 5){
  # x[1]: "original"; x[2]: "interpolated"
  falsePos = x[1] < threshold & x[2] > threshold
  falseNeg = x[1] > threshold & x[2] < threshold
  result = c(#x[1] > threshold, 
             #x[2] > threshold,
             falsePos + weightFalseNeg * falseNeg)
}
## l: location-wise, i.e. global map
absErrorMap = function(x, nout = 1){
  result = mean(abs(x[,1] - x[,2]), na.rm = TRUE)
  return(result)
}

delineationErrorMap = function(x, nout = 1, threshold = 1e-7, weightFalseNeg = 5){
  falsePos = x[,1] < threshold & x[,2] > threshold
  falseNeg = x[,1] > threshold & x[,2] < threshold
  result = mean(falsePos, na.rm = TRUE) + weightFalseNeg * mean(falseNeg, na.rm = TRUE)
  return(result)
}
