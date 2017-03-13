##########################################################################
#                     plot.Simulations                                   #
########################################################################## 
#plot.Simulations = function(x, zcol = 1, main = "", col = grey.colors(100)){
### problem: if plumes/locations have no (numeric) attributes -> error (instead it should plot without)

plot.Simulations = function(
  x, 
  ..., 
  zcol = 1, 
  main = "", 
  col = terrain.colors(100),
  maxpixels = 1e+7
){  
  
  # =================== transform values to relative =========================
  # locations properties (relative values of each numeric or factorial property)
  classesLocations = sapply(x@locations@data, class)
  classesLocations = sapply(x@locations@data, class)
  classesLocationsNum = is.element(classesLocations, c("numeric", "integer", "logical"))
  classesLocationsFac = is.element(classesLocations, c("factor"))
  for (i in which(classesLocationsFac)){
    x@locations@data[,i] = as.integer(factor(x@locations@data[,i]))
  }
  classesLocationsNF = classesLocationsNum | classesLocationsFac
  nLocationsPar = sum(classesLocationsNF) 
  if (nLocationsPar > 0){
    minLocations = lapply(x@locations@data[,classesLocationsNF, drop = FALSE], min, na.rm = TRUE)
    maxLocations = lapply(x@locations@data[,classesLocationsNF, drop = FALSE], max, na.rm = TRUE)
    propertiesLocations = (x@locations@data[,classesLocationsNF, drop = FALSE] - minLocations)/(as.list(mapply(maxLocations, minLocations, FUN = "-")))
    propertiesLocationsNA = apply(is.na(propertiesLocations), 2, all)
    propertiesLocations[propertiesLocationsNA] = 0.5
    propLocations = t(as.matrix(propertiesLocations))    
  }
  
  # plumes properties (relative values of each numeric or factorial property)
  classesPlumes = sapply(x@plumes, class)
  classesPlumesNum = is.element(classesPlumes, c("numeric", "integer", "logical"))
  classesPlumesFac = is.element(classesPlumes, c("factor"))
  for (i in which(classesPlumesFac)){
    x@plumes[,i] = as.integer(factor(x@plumes[,i]))
  }
  classesPlumesNF = classesPlumesNum | classesPlumesFac
  nPlumesPar = sum(classesPlumesNF) 
  if (nPlumesPar > 0){
    minPlumes = lapply(x@plumes[,classesPlumesNF, drop = FALSE], min, na.rm = TRUE)
    maxPlumes = lapply(x@plumes[,classesPlumesNF, drop = FALSE], max, na.rm = TRUE)
    propertiesPlumes = (x@plumes[,classesPlumesNF, drop = FALSE] - minPlumes)/(as.list(mapply(maxPlumes, minPlumes, FUN = "-")))
    propertiesPlumesNA = apply(is.na(propertiesPlumes), 2, all)
    propertiesPlumes[propertiesPlumesNA] = 0.5
    propPlumes = as.matrix(propertiesPlumes)    
  }
  
  # values (transform to matrix)
  if (canProcessInMemory(x@values)){
    propValues = matrix(getValues(subset(x@values, zcol)), nrow = ncol(x@values))    
  } else {
    propValues = subset(x@values, zcol)
  }
  
  # save the old graphics settings-- they may be needed
  def.par = par(no.readonly = TRUE)
  
  
  # =============== plotting ========================
  
  # -------------- annotations --------------------
  #par(xaxt="n", yaxt="n", bty="n", mar = c(0,0,0,0))
  par(xaxt="n", yaxt="n", bty="n", mar = c(0,0,0,0))
  # title
  par(fig = c(5/7,1,5/7,1))
  plot(0.5,0.5,type="n",ylim=c(0,1), xlim=c(0,1))
  text(0.5,0.5,paste(main), cex = 2)
  
  # annotation y-axis
  par(fig = c(0,1/7,1/7,5/7), new = TRUE)
  plot(0.5,0.5,type="n",ylim=c(0,1), xlim=c(0,1))
  text(0.5,0.5,paste("locations"), cex=1.5, srt=90)
  
  # names of locations properties
  if (nLocationsPar > 0){
    par(fig = c(5/7,1,0,1/7), new = TRUE)
    plot(0.5,0.5,type="n",ylim=c(0,1), xlim=c(0,1))
    text(labels = names(x@locations@data)[classesLocationsNF],
         x = (1:nLocationsPar)/nLocationsPar - 1/(2 * nLocationsPar), y = 1, 
         srt = 90, cex = 0.75, adj = c(1,0.5))    
  }
  
  # annotation x-axis
  par(fig = c(1/7,5/7,0,1/7), new = TRUE)
  plot(0.5,0.5,type="n",ylim=c(0,1), xlim=c(0,1))
  text(0.5,0.5,paste("plumes"), cex=1.5)
  
  # names of plumes properties
  if (nPlumesPar > 0){
    par(fig = c(0, 1/7, 5/7,1), new = TRUE)
    plot(0.5,0.5,type="n",ylim=c(0,1), xlim=c(0,1))
    text(labels = names(x@plumes)[classesPlumesNF], 
         x = 1, y = (1:nPlumesPar)/nPlumesPar - 1/(2 * nPlumesPar), 
         cex = 0.75, adj = c(1,0.5))    
  }
  
  
  #-------- data -------------------
  par(mar = c(0,0,0,0))
  
  # properties of locations
  if (nLocationsPar > 0){
    par(fig = c(5/7,1,1/7,5/7), new = TRUE)
    image(propLocations, col = col)    
  }
  
  # properties of plumes
  if (nPlumesPar > 0){
    par(fig = c(1/7,5/7,5/7,1), new = TRUE)
    image(propPlumes, col = col)    
  }
  
  # values
  par(fig = c(1/7,5/7,1/7,5/7), new = TRUE)
  
  # check if data small enough for plotting
  if (canProcessInMemory(x@values)){
    image(propValues, col = col)    
  } else {
    par(mar = c(0,0,0,0))
    par(fig = c(1/7,5/7,1/7,5/7), new = TRUE)
    image(matrix(1:4, nrow = 2), col = 0)  
    plot(propValues, col = col, add = TRUE, maxpixels = maxpixels)
         #legend = FALSE, axes = FALSE, 
         #legend.mar = 0, legend.width = 0, smallplot = c(0,0,0,0))
    warning("The values are too big to plot them all, only a random sample is plotted.")
  }
#   values = subset(x@values, zcol)
#   bs = blockSize(values, minblocks = 1)
#   if (bs$n > 1){
#     plot(propValues, col = col)
#     warning("The values are too big for this plotting method, you may plot it separately by plot(x@values).")
#   }else{
#     image(propValues, col = col) 
#   }
  box("figure", lwd = 2)
  
  # reset the graphics parameters
  par(def.par) 
}
setMethod("plot", signature(
  "Simulations"),
  definition = plot.Simulations)
# propLocations = raster(x = as.matrix(propertiesLocations), xmn = 0, xmx = 25, ymn = 0, ymx = 100)
# propPlumes = raster(t(as.matrix(propertiesPlumes)), xmn = 0, xmx = 100, ymn = 0, ymx = 25)
# propValues = raster(t(matrix(getValues(subset(simulations@values, 1)), ncol = nrow(simulations@values))), xmn = 0, xmx = 100, ymn = 0, ymx = 100)
