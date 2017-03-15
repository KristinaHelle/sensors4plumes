spplotLog = function(
  x,
  zcol,
  base,
  nTicks = 10,
#  scaleAsLog = FALSE,
  replace0 = FALSE,
  minNonzero = 0,
  ...
)
{
  # select
  if (missing(zcol)){
    zcol = names(x@data)
  }
  y = x@data[, zcol, drop = FALSE]
  
  # test, if data is numeric (no factor, logical...)
  clNum = unlist(lapply(FUN = is.numeric, X = y))
  if (any(!clNum)){
    warning(paste("Some of the data are not numeric, they cannot be transformed to logscale and are therefore excluded: ", 
                   zcol[!clNum]))
    y = y[clNum]
  }
  
  # remove negative
  y_neg = y < 0
  if (any(y_neg, na.rm = TRUE)){
    y[y_neg] = NA
    warning("The data contained negative values, these are replaced by 'NA'.")
  }
  
  # handle zero
  y0 = y == 0
  y[y0] = NA
  
  # range
  r = c(min(y, na.rm = TRUE), max(y, na.rm = TRUE))
  
  # replace values below minNonzero by 0
  if (minNonzero > r[1]){
    y_below = y < minNonzero
    y[y_below] = 0
    r[1] = 0
    warning(paste0("Some values lie between ", r[1], " and ", minNonzero, ". They are replaced by 0."))
    
    # update zero handling
    y0 = y == 0 | y0
    y[y0] = NA
    
    r = c(min(y, na.rm = TRUE), max(y, na.rm = TRUE))
  }

  
  # test bases
  if (missing(base)){
    autoBase = TRUE
    rr = r[2]/r[1]
    if (rr < 100){
      warning(paste0("The relation between maximum and minimum is ", signif(rr, 5), ":1. On a log scale there would be little difference between the values. Therefore it is plotted on normal scale."))
      base = 1
    }
    if (rr >= 100 & rr < 100000){
      base = 2
      warning(paste0("The relation between maximum and minimum is ", signif(rr,5), ":1. Data is transformed to logscale."))
    }
    if (rr >= 100000){
      base = 10
      warning(paste0("The relation between maximum and minimum is ", signif(rr,5), ":1. Data is transformed to logscale."))
    }
  } else {
    autoBase = FALSE
  }
  
  if (base != 1){
    # transform
    baseFactor = log10(base)
    yLog = log10(y)  

    
    
    # determine ticks 
    ## exponents of 10 in range
    if (autoBase & base == 2){
      baseFactor = 1
    }  
    lExp = floor(min(yLog, na.rm = TRUE)/baseFactor)
    uExp = ceiling(max(yLog, na.rm = TRUE)/baseFactor)
    tickDiff = round((uExp - lExp)/nTicks)
    
    # replace zero by next exponent below minimum
    if (replace0 & any(y0) & (class(x) != "SpatialPointsDataFrame")){
      repl0 = (lExp - 2 * tickDiff) * baseFactor
      if (autoBase & base == 2){
        repl0 = lExp - 2 * log10(2)
      }
      yLog[y0] = repl0 
    } 
    if (!any(y0)){
      replace0 = FALSE
    }
    
    if (autoBase & base == 2){
      AtAll = c(1,2,5) * 10^rep((lExp - 1):(uExp + 1), each = 3)
    } else {
      if (tickDiff == 0){
        tickDiff = signif((uExp - lExp)/nTicks, 2)
        warning("The chosen 'base' and 'nTicks' do not result in nice numbers.")
      }
      AtAll = base^seq(lExp - tickDiff, uExp + tickDiff, tickDiff)
    }

    AtIn = range(which(AtAll >= min(y, na.rm = TRUE) & AtAll <= max(y, na.rm = TRUE)))
    
    Breaks = log10(AtAll[(AtIn[1] - 1):(AtIn[2] + 1)])
    nBr = length(Breaks)
    Labels = AtAll[AtIn[1]:AtIn[2]]
    At = log10(Labels)
    
    # correct labelling of (replaced) zero
    if (replace0 & (class(x) != "SpatialPointsDataFrame")){
      Breaks = c(repl0, Breaks)
      At = c(repl0, At)
      Labels = c(0, Labels)
      nBr = nBr + 1
    } 
    
    # spatial object with log scaled values
    z = x
    z@data = yLog

    
    # plot
    if (class(x) != "SpatialPointsDataFrame"){
      if (is(x, "SpatialIndexDataFrame") & length(z@data) == 1){ # extra treatment needed, else no plotting
        zDataNames = names(z@data)
        zData = cbind(z@data, Q = 1)
        z@data = zData
        spplot(z,
               zcol = zDataNames,
               colorkey = list(
                 labels=list(
                   at = At,
                   labels = Labels
                 )), 
               at = Breaks,
               ...) 
      } else {
        spplot(z,
               colorkey = list(
                 labels=list(
                   at = At,
                   labels = Labels
                 )), 
               at = Breaks,
               ...)  
      }
    } else {
      # SPointsDF: plot original data with changed cuts
      Cuts = 10^Breaks
      if (replace0 & any(y0)){
        y[y0] = 0
        Cuts = c(0, Cuts)
      } 
      z = x
      z@data = y
 
      spplot(z,
             cuts = Cuts,
             ...)
    }
  } else {
    if (replace0 & any(y0)){
      y[y0] = 0
    } 
    z = x
    z@data = y
    if (is(x, "SpatialIndexDataFrame") & length(z@data) == 1){ # extra treatment needed, else no plotting
      zDataNames = names(z@data)
      zData = cbind(z@data, Q = 1)
      z@data = zData
      spplot(z,
             zcol = zDataNames,
             ...) 
    } else {
      spplot(z, ...) 
    }
  }
}











# spplotLog = function(
#   x,
#   zcol,
#   base = 10,
#   nTicks = 10,
#   scaleAsLog = FALSE,
#   replace0 = TRUE,
#   minNonzero = 0,
#   ...
# ){
#   if (missing(zcol)){
#     zcol = names(x@data)
#   }
#   y_0 = x@data[, zcol, drop = FALSE]
#   y_0[y_0 <= minNonzero] = 0
#   # transform to logscale
#   y = log(y_0)/log(base)
#   if (replace0){
#     # replace zero by next exponent below minimum
#     y[y_0 == 0] = NA
#     minY = min(y, na.rm = TRUE)  
#     y[y_0 == 0] = floor(minY - 1)
#   }
#   
#   # find nice levels for color breaks and ticks (try to reproduce what spplot does)
#     if (all(is.finite(na.omit(unlist(y))))){ # range limited? else it cannot be split into equal parts
#       r = diff(range(y, na.rm = TRUE))
#     } else {
#       spplot(x)
#       warning("Plot on original scale (not logarithmic) because it contains values < = 0 and 'replace0 = FALSE' or other reason to make tranformed values unlimited.")
#     } 
#   
#     r_B = r/nTicks # about desired tick distance (for nTicks ticks)
#     rExp = floor(log(r_B)/log(10))
#     rMan = ceiling(r_B/(10^rExp))
#     ticksDiff = rMan * 10^rExp # desired tick distance: "nice" value
#     
#     At = seq(floor(min(y, na.rm = TRUE)/ticksDiff) * ticksDiff,
#              floor(max(y, na.rm = TRUE)/ticksDiff) * ticksDiff,
#              ticksDiff) 
#     if (!scaleAsLog){
#       Labels = signif(base^At, 3)
#     } else {
#       Labels = paste0(base, "^", At)  
#     }
#     if (replace0 & any(y_0 == 0, na.rm = TRUE)){ # correct labelling of (replaced) zero
#       subzeroTicks = floor(minY - 1) >= At
#       Labels = Labels[!subzeroTicks]
#       At = At[!subzeroTicks]
#       Labels = c("0", Labels)
#       At = c(floor(minY - 1), At) 
#     }
#     
#     # spatial object with log scaled values
#     z = x
#     z@data = y
#     
#     # plot
#     spplot(z,
#            colorkey = list(
#              labels=list(
#                at = At,
#                labels = Labels
#              )), ...)
# } 



#paste0(base, "^", At)  
# else {
# 
# rLog = diff(range(yLog, na.rm = TRUE))
# r_B = rLog/nTicks # about desired tick distance (for nTicks ticks)
# rExp = floor(log(r_B)/log(10))
# rMan = ceiling(r_B/(10^rExp))
# ticksDiff = rMan * 10^rExp # desired tick distance: "nice" value
# 
# 
# 
#   
#   # break color at these values
#   At = seq(floor(min(yLog, na.rm = TRUE)/ticksDiff) * ticksDiff,
#            ceiling(max(yLog, na.rm = TRUE)/ticksDiff) * ticksDiff,
#            ticksDiff) 
#   
#   # labels  
#   if (!scaleAsLog){
#     Labels = signif(base^At, 3)
#   } else {
#     Labels = expression(paste0(base, "^", At))  
#   }      
# }

# correct labelling of (replaced) zero
# if (replace0){ 
#   subzeroTicks = floor(minY - 1) 
#   Labels = Labels[!subzeroTicks]
#   At = At[!subzeroTicks]
#   Labels = c("0", Labels)
#   At = c(floor(minY - 1), At) 
# }


    
