##################################################
#            changeSimulationsPath               #
#               copySimulations                  #
##################################################

changeSimulationsPath = function(simulations, path){
  if (is(simulations, "Simulations")){
    raster = simulations@values
  } else {
    raster = simulations
  }  
  
  if (!inMemory(raster)){
    if (length(path) == 1){
      if (is(raster, "RasterLayer")){
        newRaster = raster(path)
      } else {
        newRaster = brick(path)
      }
    } else {
      newRaster = stack(path)
    } 
    if (is(simulations, "Simulations")){
      newSimulations = Simulations(
        locations = simulations@locations,
        plumes = simulations@plumes,
        values = newRaster)
    } else {
      newSimulations = newRaster
    }  
  } else {
    newSimulations = simulations
  }
  return(newSimulations)
}
  

copySimulations = function(simulations, 
                           newPath, 
                           #newName, 
                           newFile, 
                           overwrite = FALSE,
                           deleteOld = FALSE
                           ){
  if (!inMemory(simulations@values) & !missing(newPath)){
    if (!file.exists(newPath)){
      stop(paste0("There is no directory '", newPath, "'."))
    }
    # determine files
    nLay = nlayers(simulations@values)
    oldPaths = character(nLay)
    for (i in 1:nLay){
      oldPaths[i] = simulations@values[[i]]@file@name
    }
    layerNames = sapply(lapply(X = strsplit(oldPaths, "/"), FUN = rev), "[[", 1)
    layerNames = strsplit(layerNames, "[.]")
    suffixes = sapply(layerNames, "[[", 2)
    fileNames =  sapply(layerNames, "[[", 1)
    
    # copy files
    for (i in 1:nLay){
      file.copy(from = oldPaths[i], to = newPath, overwrite = overwrite) # .grd files
      if (suffixes[i] == "grd"){
        file.copy(from = paste0(strsplit(oldPaths, ".grd")[[i]], ".gri"), to = newPath, overwrite = overwrite) # .gri files        
      }
    }
    
    # update paths
    simulationsNew = changeSimulationsPath(simulations,
                                           path = paste0(newPath, "/", paste0(fileNames, ".", suffixes)))
    simulations = simulationsNew
    
    #for (i in 1:nLay){
    #  simulations@values[[i]]@file@name = paste0(newPath, "/", layerNames[i], ".grd")
    #}
    
    # delete old files
    if (deleteOld){
      file.remove(oldPaths)
      if (suffixes[i] == ".grd"){
        file.remove(paste(strsplit(oldPaths, ".grd"), ".gri", sep = ""))
      }  
    }
  }  
    
  if (!missing(newFile)){ # save
    if (!missing(newPath)){
      save(simulations, file = paste0(newPath, "/", newFile, ".Rdata"))
    } else {  
      save(simulations, file = paste0(newFile, ".Rdata"))
    }
  }
  return(simulations)
}

