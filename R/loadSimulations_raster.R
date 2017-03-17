loadSimulations_raster = function(
  basicPath,
  filePaths,
  region,
  bBox,
  nameSave = NA,
  overwrite = FALSE,
  nP,
  nL,
  nK,
  multilayer,
  ...)
{
  print("Load data...")
  
  simulations = list()
  if (nK == 1){# 1 kind
    if (is.na(multilayer)){# use only 1st layer
      simulations[[1]] = stack(paste0(basicPath, "/", filePaths[,1]), bands = 1) 
    } else {
      if (multilayer == "plumes"){# use all layers
        simulations[[1]] = stack(paste0(basicPath, "/", filePaths[,1]))   
        nP = nlayers(simulations[[1]])
      }   
      if (multilayer == "kinds"){
        stop("ERROR nK == 1 & multilayer == 'kinds'")
      }      
    }
  } else {
    if (is.na(multilayer)){# kinds in different files
      paths5 = split(filePaths[,1], filePaths[,3])
      for (i in seq(along = paths5)){
        simulations[[names(paths5)[i]]] = stack(paste0(basicPath, "/", paths5[[i]]))
      }
    } else {
      if (multilayer == "kinds"){
        simulations0 = stack(paste0(basicPath, "/", filePaths[,1])) 
        if (nlayers(simulations0) != nK * nP){
          stop(paste0("Wrong number of layers in some file; there should be ", nP * nK, " layers in all files together, but there are ", nlayers(simulations0),"."))
        }
        for (i in 1:nK){
          simulations[[i]] = subset(simulations0, subset = seq(i, nK * nP, nK))
        }
        if (strsplit(filePaths[1,1], "\\.")[[1]][2] == "grd"){
          simulations00 = brick(paste0(basicPath, "/", filePaths[1,1]))
          names(simulations) = names(simulations00)
        } else {
          names(simulations) = 1:nK
        }
        
      }
      if (multilayer == "plumes"){
        stop("ERROR nK > 1 & multilayer == 'plumes'")
      }      
    }
    
  }
  print("All data loaded.")
  
  # if bBox is given, cut out this part
  if(!missing(bBox)){
    for (i in 1:nK){
      cropped = crop(simulations[[i]], bBox)
      simulations[[i]] = cropped
    }   
  } 
  #  for (i in 1:nK){
  #    writeRaster(simulations[[i]], filename = paste(nameSave, "/cropped_", i, ".grd", sep = ""),
  #                overwrite = overwrite)
  #  }
  
  values = list()
  for (i in 1:nK){
    values[[i]] = raster(nrows = nL, ncols = nP,
                         xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                         crs = "+init=epsg:4326")
  }
  
  names(values) = names(simulations)[1:nK]
  
  # transform 'simulations' to matrices with row ~ location, column ~ plume
  bs = blockSize(simulations[[1]], minblocks = 1, n = nP) 
  print(paste("Data is split into ", bs$n, " block(s) for processing.", sep = ""))
  for (i in 1:nK){
    dataType_i = dataType(simulations[[i]])
    values[[i]] = writeStart(values[[i]], 
                             filename = paste(nameSave, "_", i, ".grd", sep = ""), 
                             overwrite = overwrite,
                             dataType = dataType_i)
    k = 1
    rowlength = dim(simulations[[i]])[2]
    for (j in 1:bs$n){ 
      in_j = getValues(simulations[[i]], row = bs$row[j], nrows = bs$nrows[j])     
      if(!is.null(region)){
        locations_j = (bs$row[j] - 1) * rowlength  + 1:(bs$nrows[j] * rowlength)
        select_j = region[locations_j]
        in_j = in_j[select_j,, drop = FALSE]
        k_j = sum(select_j)
      }else{
        k_j = dim(in_j)[1]
      }
      values[[i]] = writeValues(x = values[[i]], 
                                v = t(in_j), 
                                start = k) 
      k = k + k_j
      print(paste("Block ", j, " processed for kind ", names(values)[i], "."))
    }
  }  
  
  for (i in 1:nK){
    values[[i]] = writeStop(values[[i]])
  }
  # if crop has been saved, delete these intermediate files
  #  if (!missing(bBox)){
  #    file.remove(paste(cropped[[2]], ".grd", sep = ""))
  #    file.remove(paste(cropped[[2]], ".gri", sep = ""))            
  #  }
  # turn into one Raster* object, kinds ~ layers
  Values = stack(values)
  print("All values processed.")
  
  # extract locations
  Locations = as.data.frame(coordinates(simulations[[1]]))
  if(!is.null(region)){
    Locations = Locations[region,, drop = FALSE]
  }  
  Locations$index = 1:nL
  coordinates(Locations) = 1:2
  proj4string(Locations) = CRS(proj4string(simulations[[1]]))
  gridded(Locations) = TRUE
  #  Locations = as(Locations, "SpatialPolygridDataFrame")
  
  # extract plume properties
  if (is.na(multilayer)){
    if (dim(filePaths)[2] > 2){
      Plumes = data.frame(name = split(filePaths[,2], filePaths[,3])[[1]])      
    } else {
      Plumes = data.frame(name = filePaths[,2]) 
    }     
  } else {
    if (multilayer == "plumes"){
      Plumes = data.frame(name = names(simulations[[1]]))       
    } else {
      Plumes = data.frame(name = filePaths[,2])       
    }
  }
  
  
  # turn into Simulations
  out = Simulations(locations = Locations, 
                    plumes = Plumes,
                    values = Values)
  
  return(out)  
}