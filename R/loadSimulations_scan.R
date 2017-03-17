loadSimulations_scan = function(
  basicPath,   
  filePaths,
  region,
  nameSave = NA, 
  overwrite,
  nP,
  nL,
  nK,
  multilayer,
  ...)
{
  print("Load data...")  
  
  # make raster objects to save data: row ~ plume, col ~ location  
  if (!is.na(multilayer)){
    if (multilayer == "plumes"){
      data0 = data.frame(x = 0)
      for (i in 1:nrow(filePaths)){
        data0 = cbind(data0, read.table(file = paste0(basicPath, "/", filePaths[i,1]), 
                                        header = TRUE, nrows = 1, ...))
      } 
      data0[,1] = NULL
      nP = ncol(data0) 
    } else {
      data0 = read.table(file = paste0(basicPath, "/", filePaths[1,1]), 
                         header = TRUE, nrows = 1, ...)
    }
  }

  simulations = list()
  for (i in 1:nK){    
    simulations[[i]] = raster(nrows = nP, ncols = nL,
                              xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                              crs = "+init=epsg:4326")
    simulations[[i]] = writeStart(simulations[[i]],
                                  filename = paste(nameSave, "_switched_", i, ".grd", sep = ""), 
                                  overwrite = overwrite)
                                  #, dataType = dataType(data1[,1]))
  }
  kindCount = rep(0, nK)
  if (dim(filePaths)[2] > 2){
    fileKind = as.integer(factor(filePaths[,3]))
  }

  
  for (j in 1:nrow(filePaths)){
    data1 = scan(file = paste0(basicPath, "/", filePaths[j,1]), ...)
    
    if (!is.null(region)){
      nLin = length(region)
      if (length(data1) %% nLin != 0){
        stop(paste0("The numbers of values in the data files differ, such data cannot be used. The error occured for: ", 
                    filePaths[j,1]))
      }
      if (!is.na(multilayer)){
        if (multilayer == "kinds"){
          if (length(data1) != nLin * nK){
            stop(paste0("The numbers of values in file ", filePaths[j,1], " was not ", nLin * nK, "."))           
          }
        }
      }
      data2 = matrix(data1, byrow = TRUE, nrow = nLin) 
      data = data2[region,, drop = FALSE]
    } else {
      if (length(data1)%%nL != 0){
        stop(paste0("The numbers of values in the data files differ, such data cannot be used. The error occured for: ", 
                    filePaths[j,1]))
      }      
      if (!is.na(multilayer)){
        if (multilayer == "kinds"){
          if (length(data1) != nL * nK){
            stop(paste0("The numbers of values in file ", filePaths[j,1], " was not ", nL * nK, "."))            
          }
        }
      }
      
      data = matrix(data1, byrow = TRUE, nrow = nL) 
    }
    
    if (nK == 1){# 1 kind
      if (is.na(multilayer)){# use only 1st layer
        simulations[[1]] = writeValues(x = simulations[[1]], 
                                       v = data[,1],
                                       start = j)
      } else {
        if (multilayer == "plumes"){# use all layers
          for (k in 1:ncol(data)){
            simulations[[1]] = writeValues(x = simulations[[1]], 
                                           v = data[,k],
                                           start = kindCount[1] + k)
          }
          kindCount[1] = kindCount[1] + ncol(data)    
        }   
        if (multilayer == "kinds"){
          stop("ERROR nK == 1 & multilayer == 'kinds'")
        }      
      }
    } else {
      if (is.na(multilayer)){# kinds in different files
        kindCount[fileKind[j]] = kindCount[fileKind[j]] + 1
        simulations[[fileKind[j]]] = writeValues(x = simulations[[fileKind[j]]], 
                                                      v = data[,1],
                                                      start = kindCount[fileKind[j]])
      } else {
        if (multilayer == "kinds"){
          for (i in 1:nK){
            kindCount[i] = kindCount[i] + 1
            simulations[[i]] = writeValues(x = simulations[[i]], 
                                           v = data[,i],
                                           start = kindCount[i])
          }
        }
        if (multilayer == "plumes"){
          stop("ERROR nK > 1 & multilayer == 'plumes'")
        }      
      }  
    }
  }
  
  for (i in 1:nK){
    simulations[[i]] = writeStop(simulations[[i]])
  }
  names(simulations) = 1:nK
  if (!is.na(multilayer)){
    if (multilayer == "kinds"){
      names(simulations) = names(data0)      
    }
  }
  if (dim(filePaths)[2] > 2){
    names(simulations) = sort(unique(filePaths[,3]))
  }  
  print("All data loaded.")
  

  # transpose (writes to temporary file first, it rather should write to the given file directly)
  values = list()
  for (i in 1:nK){
    values[[i]] = t(simulations[[i]])         
#   values[[i]] = writeRaster(values[[i]],
#                              filename = paste(nameSave, "simulations_", i, ".grd", sep = ""), 
#                              overwrite = overwrite)
  }
  names(values) = names(simulations)
  Values = stack(values)
  
  # extract plume properties
  if (is.na(multilayer)){
    if (dim(filePaths)[2] > 2){
      Plumes = data.frame(name = split(filePaths[,2], filePaths[,3])[[1]])      
    } else {
      Plumes = data.frame(name = filePaths[,2]) 
    }     
  } else {
    if (multilayer == "plumes"){
      Plumes = data.frame(name = names(data0))       
    } else {
      Plumes = data.frame(name = filePaths[,2])       
    }
  }
  
  
  # remove intermediate files
  rm(simulations)
  for (i in 1:nK){
    file.remove(paste(nameSave, "_switched_", i, ".grd", sep = ""))  
    file.remove(paste(nameSave, "_switched_", i, ".gri", sep = ""))    
  }
  
  # return as list
  out = list(plumes = Plumes,
             values = Values)
  
  return(out) 
}
