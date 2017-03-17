######################################################################
# cbind.Simulations (by plume)                                       #
######################################################################

cbind.Simulations = function(
  ..., # Simulations objects; if locations are the same and as well names of values and plumes, plumes (and values) are combined 
  #  saveName, # filename for new file to be created, if missing, random unique name is used
  #  saveDir, 
  nameSave = NA, 
  overwrite = FALSE
){
  dots = list(...)
  n = length(dots)
  if (n >= 1){   
    nL = nrow(dots[[1]]@values)
    if (n == 1){
      return(dots[[1]])
    }
  }else{
    stop("No object(s) given! Returns NULL.")
    newSim = NULL
  }
  
  # ----------- test if everything is identical ------------------ #
  if (n >= 2){
    for (l in 2:n){
      if (!identical(dots[[1]]@locations, dots[[l]]@locations)){
        stop(paste("Simulations cannot be combined as their 'locations' differ. The",l, "th object differs from the previous ones."))
      }
      layerNames = names(dots[[1]]@values)
      if (!identical(layerNames, names(dots[[l]]@values))){
        stop(paste("Simulations cannot be combined as the layers of their 'values' differ by name. The",l, "th object differs from the previous ones."))
      }
      nLay = length(layerNames)
      if (!identical(names(dots[[1]]@plumes), names(dots[[l]]@plumes))){
        stop(paste("Simulations cannot be combined as their 'plumes' differ by name. The",l, "th object differs from the previous ones."))
        
      }
    }  
  }
  
  # if no saveName given, create random name extension
  #  if (missing(saveDir)){
  #    saveDir = getwd()
  #    warning(paste("No 'saveDir' indicated, files are saved in current directory:", getwd()))
  #  }
  #  if (missing(saveName)){
  #    randName = tempfile(pattern = paste("simulations_", layerNames[1], "_", sep = ""), tmpdir = saveDir, fileext = ".grd")
  #    saveName = strsplit(strsplit(randName, paste(getwd(), "\\\\", "simulations_", layerNames[1], "_", sep = ""))[[1]][2], ".grd")[[1]]
  #  }
  
  # number of plumes per object
  nP = integer(n)
  for (l in 1:n){
    nP[l] = ncol(dots[[l]]@values)
  }
  
  #----------- combine values ---------------#
  bs = blockSize(raster(nrow = nL, ncol = sum(nP)), minblocks = 1, n = n)
  values_new = list()
  for (i in seq(along = layerNames)){
    values_new[[layerNames[i]]] = raster(nrow = nL, ncol = sum(nP), xmn = -90, xmx = 90, ymn = -90, ymx = 90, crs = "+init=epsg:4326")
    if(bs$n > 1){
      values_new[[layerNames[i]]] = writeStart(values_new[[layerNames[i]]], 
                                               filename = paste0(nameSave, "_", layerNames[i], ".grd"),
                                               overwrite = overwrite)
      
      #        filename = paste(saveDir, "/simulations_", layerNames[i], "_", saveName, ".grd", sep = ""), overwrite = overwrite)
      #       overwritten = "overwritten"
      #       if (!overwrite){
      #         overwritten = "not overwritten"
      #       }
      #       warning(paste0("The resulting raster is saved as ",
      #                      paste0(nameSave, "_", layerNames[i], ".grd"),
      #                      "existing files are ", overwritten, "."))
    }
  }  
  cnP = cumsum(c(0,nP))
  if (bs$n <= 1){
    in_j = matrix(nrow = nL * sum(nP), ncol = nLay, byrow = FALSE) 
    for (l in 1:n){
      in_jl = getValues(dots[[l]]@values) 
      index = rep(cnP[l] + sum(nP) * (0:(nL - 1)), each = nP[l]) + rep(1:nP[l], times = nL)
      in_j[index,] = in_jl
    }
    for (i in seq(along = layerNames)){
      values(values_new[[layerNames[i]]]) = in_j[,i]
    }
  }else{
    for (j in 1:bs$n){
      # get and order values
      in_j = matrix(nrow = bs$nrows[j] * sum(nP), ncol = nLay, byrow = FALSE)
      for (l in 1:n){
        in_jl = getValues(dots[[l]]@values, row = bs$row[j], nrows = bs$nrows[j]) 
        index = rep(cnP[l] + sum(nP) * (0:(bs$nrows[j]-1)), each = nP[l]) + rep(1:nP[l], times = bs$nrows[j])
        in_j[index,] = in_jl
      }
      # write values to file
      for (i in seq(along = layerNames)){
        values_new[[layerNames[i]]] = writeValues(x = values_new[[layerNames[i]]], v = in_j[,i], start = bs$row[j])
      }
    }  
    for (i in seq(along = layerNames)){
      values_new[[layerNames[i]]] = writeStop(values_new[[layerNames[i]]])
    }     
  }               
  
  Values_new = stack(values_new)
  
  # --------------- combine plumes ----------- #
  Plumes_new = dots[[1]]@plumes
  for (l in 2:n){
    Plumes_new = rbind(Plumes_new, dots[[l]]@plumes)
  }
  
  # output
  newSim = Simulations(
    locations = dots[[1]]@locations,
    plumes = Plumes_new,
    values = Values_new)
  return(newSim)
}



