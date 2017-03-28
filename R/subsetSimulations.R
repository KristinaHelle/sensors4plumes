#############################################################
#                   subset.Simulations                      #
#############################################################
# selectRow # internal
# selectCol # internal

# subset

# helper functions
selectCol = function(colSelect, nrow, ncol){
  colSelectIn = colSelect[colSelect >0 & colSelect <= ncol]
  if (length(colSelectIn) < length(unique(colSelect))){
    warning("Some of the chosen indices are outside, they are ignored.")
  }
  return(rep(colSelectIn, times = nrow) + ncol * rep(0:(nrow-1), each = length(colSelectIn)))
}
selectRow = function(rowSelect, nrow, ncol){
  rowSelectIn = rowSelect[rowSelect > 0 & rowSelect <= nrow]
  if (length(rowSelectIn) < length(unique(rowSelect))){
    warning("Some of the chosen indices are outside, they are ignored.")
  }
  return(ncol * rep(rowSelectIn - 1, each = ncol) +  rep(1:ncol, times = length(rowSelectIn)))  
}


subset.Simulations = function(
  x, # Simulations object
  ..., # for consistency with generic 'subset'
  locations , #= 1:length(simulations@locations) index of simulations@locations to extract
  plumes , # 1:nrow(simulations@plumes)index of simulations@plumes to extract
  kinds , # 1:nlayers(simulations@values)index of simulations@values to extract
  dataLocations , #  1:ncol(simulations@locations@data)columns of simulations@locations@data
  dataPlumes, #1:ncol(simulations@plumes)columns of simulations@plumes
  nameSave = NA,
#  saveName, # where to save the file if necessary
#  saveDir = ".",
  overwrite = FALSE,
  valuesOnly = FALSE # then only the values are subsetted and returned; only in this case multiple locations and plumes are taken into account; else they are ignored
){
  # - - - - subsetting that does not require splitting of the values - - - -
  # subset values by kind
  if(!missing(kinds)){  
    x@values = subset(x@values, kinds, drop = FALSE)
  }
  # properties of (remaining) values
  layerNames = names(x@values)
  nL0 = nrow(x@values)
  nP0 = ncol(x@values)

  # subset columns of x@locations@data and x@plumes (not for 'valuesOnly')  
  if (!valuesOnly){
    if(!missing(dataLocations)){
      x@locations@data = x@locations@data[,dataLocations, drop = FALSE]
    }
    if(!missing(dataPlumes)){
      x@plumes = x@plumes[,dataPlumes, drop = FALSE]
    }      
  }  

  # - - - -  subsetting of values by locations and plumes - - - - - 
  data = x@values
  # correct locations 
  isLoc = FALSE
  if(!missing(locations)){   
    locations[locations <= 0 | locations > nL0] = NA
    if (!valuesOnly){
      locations = unique(locations[!is.na(locations)])
    } else {
      locationsNA = is.na(locations)
      locations = locations[!locationsNA]
      nLNA = length(locationsNA)
    }
    nL = length(locations)
    if(!identical(locations, 1:nL0)){      
      isLoc = TRUE
    }
  }else{
    nL = nL0
  }
  if(nL == 0){
    stop("There are 0 valid 'locations' selected.")
  }
  
  # correct plumes
  isPl = FALSE
  if (!missing(plumes)){
    plumes[plumes <= 0 | plumes > nP0] = NA
    if (!valuesOnly){
      plumes = unique(plumes[!is.na(plumes)])
    } else {
      plumesNA = is.na(plumes)
      plumes = plumes[!plumesNA]
      nPNA = length(plumesNA)
    }
    nP = length(plumes)
    if(!identical(plumes, 1:ncol(x@values))){
      isPl = TRUE
    }
  }else{
    nP = nP0
  }
  if(nP == 0){
    stop("There are 0 valid 'plumes' selected.")
  }    
    
  #if (class(x@values) != "RasterStack"){
  #  stop("The x@values have to be of class 'RasterStack'.")
  #}
  
  # subset plumes and locations
  ## subset values
  if (isLoc || isPl){
    # input: needs processing in chunks?   
    bs = blockSize(data, minblocks = 1, n = 2)
    # results: needs to be saved to disk?
    if (!valuesOnly){
      bs_out = blockSize(raster(nrow = nL, ncol = nP), n = length(layerNames) * 2, minblocks = 1)      
    } else {
      bs_out = blockSize(raster(nrow = nLNA, ncol = nPNA), n = length(layerNames) * 2, minblocks = 1)
    }

    if(bs_out$n > 1){# then saving to disk is necessary
      # if no saveName given, create random name extension
      if (is.na(nameSave)){
        stop("Data needs to be saved to disk as they are too big to be processed in memory, please indicate 'nameSave'.")
        # randName = tempfile(pattern = paste("simulations_", layerNames[1], "_", sep = ""), tmpdir = saveDir, fileext = ".grd")
        # saveName = strsplit(strsplit(randName, paste(getwd(), "\\\\", "simulations_", layerNames[1], "_", sep = ""))[[1]][2], ".grd")[[1]]
      }  
      warning(paste("Values are saved to disk:", nameSave, "_", layerNames, ".grd; ", sep = ""))
    }
    
    values_new = list()
    for (i in seq(along = layerNames)){
      if (!valuesOnly){
        values_new[[layerNames[i]]] = raster(nrow = nL, ncol = nP, 
                                             xmn = -90, xmx = 90, ymn = -90, ymx = 90, 
                                             crs = "+init=epsg:4326")        
      } else {
        values_new[[layerNames[i]]] = raster(nrow = nLNA, ncol = nPNA, 
                                             xmn = -90, xmx = 90, ymn = -90, ymx = 90, 
                                             crs = "+init=epsg:4326")  
      }
      if(bs$n <= 1){# input loaded at once
        in_j = getValues(data[[layerNames[i]]])     
        nR = nL0
        nC = nP0
        if (isLoc){        
          in_j = in_j[selectRow(locations, nR, nC)]
          nR = nL
        }
        if (isPl){
          in_j = in_j[selectCol(plumes, nR, nC)]
        }        
        if (valuesOnly){
          In_j = rep(1, nLNA * nPNA)
          if (isLoc){
            In_j[selectRow(which(locationsNA), nLNA, nPNA)] = NA
          }
          if (isPl){
            In_j[selectCol(which(plumesNA), nLNA, nPNA)] = NA
          }
          In_j[!is.na(In_j)] = in_j
          in_j = In_j
        }
        values(values_new[[layerNames[i]]]) = in_j
      }else{# input read in chunks
        if (bs_out$n > 1){ # in this case create new file for output
          values_new[[layerNames[i]]] = 
            writeStart(values_new[[layerNames[i]]],
                       filename = paste0(nameSave, "_", layerNames[i], ".grd"), overwrite = overwrite)          
        }
        k = 1
        for (j in 1:bs$n){
          in_j = getValues(x@values[[layerNames[i]]], row = bs$row[j], nrows = bs$nrows[j])     
          nR = bs$nrows[j]
          nC = nP0
          indicesTheseLoc = bs$row[j] - 1 + 1:bs$nrows[j]
          if (isLoc){  
            locInChunk = (locations >= bs$row[j]) & (locations < (bs$row[j] + bs$nrows[j]))
            theseLoc = locations[locInChunk] - (bs$row[j] - 1) 
            indicesTheseLoc = (1:nL)[locInChunk]         
            in_j = in_j[selectRow(theseLoc, nR, nC)]
            nR = length(theseLoc)
          } else {
            theseLoc = 1:bs$nrows[j] 
            nR = length(theseLoc)
          }
          if (nR > 0){
            if (isPl){
              in_j = in_j[selectCol(plumes, nR, nC)]
            }
            if(bs_out$n > 1){
              if (valuesOnly){
                nLThese = length(indicesTheseLoc)
                In_j = rep(1, nLThese * nPNA)
                if (isPl){
                  In_j[selectCol(which(plumesNA), nLThese, nPNA)] = NA
                }
                In_j[!is.na(In_j)] = in_j
                in_j = In_j
              }  
              for (l in seq(along = theseLoc)){
                if (!valuesOnly){
                  values_new[[layerNames[i]]] = writeValues(x = values_new[[layerNames[i]]], 
                                                            v = in_j[(l-1) * nP + 1:nP], 
                                                            start = indicesTheseLoc[l])                   
                } else {
                  values_new[[layerNames[i]]] = writeValues(x = values_new[[layerNames[i]]], 
                                                            v = in_j[(l-1) * nPNA + 1:nPNA], 
                                                            start = which(!locationsNA)[indicesTheseLoc[l]]) 
                }       
              }
              if (valuesOnly){
                for (l in 1:sum(locationsNA)){
                  values_new[[layerNames[i]]] = writeValues(x = values_new[[layerNames[i]]], 
                                                            v = rep(NA, nPNA),
                                                            start = which(locationsNA)[l]) 
                }                
              }
            }else{
              if (valuesOnly){
                nLThese = length(indicesTheseLoc)
                In_j = rep(1, nLThese * nPNA)
                if (isPl){
                  In_j[selectCol(which(plumesNA), nLThese, nPNA)] = NA
                }
                In_j[!is.na(In_j)] = in_j
                in_j = In_j
                values_new[[layerNames[i]]][which(!locationsNA)[indicesTheseLoc],] = in_j 
              } else {
                values_new[[layerNames[i]]][indicesTheseLoc,] = in_j  
              }
            }
            k = k + nR          
          }
        }  
        if(bs_out$n > 1){
          values_new[[layerNames[i]]] = writeStop(values_new[[layerNames[i]]])          
        }
      }  
    }  
    Values_new = stack(values_new)
    
    # put together output
    if (!valuesOnly){
      ## subset locations
      if (isLoc){
        locations_new = subsetSDF(x = x@locations, locations = locations) 
      }
      
      ## subset plumes
      if (isPl){
        plumes_new = x@plumes[plumes,, drop = FALSE]  
      }
      
      # combine to new x object
      if (isLoc && isPl){
        simulations_new = Simulations(
          locations = locations_new,
          plumes = plumes_new,
          values = Values_new)
      }
      
      if (isLoc && !isPl){
        simulations_new = Simulations(
          locations = locations_new,
          plumes = x@plumes,
          values = Values_new)      
      }
      if (!isLoc && isPl){
        simulations_new = Simulations(
          locations = x@locations,
          plumes = plumes_new,
          values = Values_new)
      }
      
      if (!isLoc && !isPl){
        simulations_new = x
      }      
    } else {
      simulations_new = Values_new
    }
  } else {
    if (!valuesOnly){
      simulations_new = x
    } else {
      simulations_new = x@values
    }
  }
  return(simulations_new)
}
setMethod("subset", signature("Simulations"), definition = subset.Simulations)
