###########################################################
# postprocessing                                          #
###########################################################
SDLonLat = function(simulations, SD){
  # turn SD into list
  if (!is.list(SD)){
    if (is.matrix(SD)){
      SDs = list()
      for (i in 1:nrow(SD)){
        SDs[[i]] = SD[i,]
      }
      SD = SDs
    } else {# single SD
      SD = list(SD)      
    }
  }
  coordOptLonLat = list()
  for (i in seq(along = SD)){
    coordOpt = as.data.frame(coordinates(simulations@locations)[SD[[i]],])
    coordinates(coordOpt) = 1:2
    proj4string(coordOpt) = CRS(proj4string(simulations@locations)) 
    coordOptLonLat[[i]] = spTransform(coordOpt, CRS("+init=epsg:4326"))
  }
  return(coordOptLonLat)
}

# - -- - - - - - - - compare SDs  - - - - - - - - - - -- 
# spatial similarity
similaritySD = function(simulations, SD, referenceSD, type = "equal", k = 9, ...){
  # turn SD into list
  if (!is.list(SD)){
    if (is.matrix(SD)){
      SDs = list()
      for (i in 1:nrow(SD)){
        SDs[[i]] = SD[i,]
      }
      SD = SDs
    } else {# single SD
      SD = list(SD)      
    }
  }
  if (type == "EarthMoversDistance"){
    # generate matrix that represents SD
    if (! is.element(class(simulations@locations),
                     c("SpatialPixelsDataFrame", 
                       "SpatialPolygridDataFrame"))){
      stop("'EarthMoversDistance' can only be used with 'simulations' that have 'locations' that are 'Pixels' or 'Polygrid'.")
    } else {
      switch(class(simulations@locations),
             "SpatialPixelsDataFrame" = {
               simulations@locations@data$index = 1:nLocations(simulations)
               simGrid = as(simulations@locations, "SpatialGridDataFrame")
               simMatrix = matrix(simGrid@data$index, byrow = TRUE,
                                  nrow = simulations@locations@grid@cells.dim[2],
                                  ncol = simulations@locations@grid@cells.dim[1])
               referenceMatrix = simMatrix
               referenceMatrix[!is.element(referenceMatrix, referenceSD) & !is.na(referenceMatrix)] = 0
               referenceMatrix[is.element(referenceMatrix, referenceSD)] = 1
             },
             "SpatialPolygridDataFrame" = {
               simMatrix = matrix(simulations@locations@index, byrow = TRUE,
                                  nrow = simulations@locations@grid@cells.dim[2],
                                  ncol = simulations@locations@grid@cells.dim[1])
               referenceMatrix = simMatrix
               referenceArea = as.data.frame(table(simMatrix[is.element(simMatrix, referenceSD)]))
               referenceMatrix[!is.element(referenceMatrix, referenceSD) & !is.na(referenceMatrix)] = 0
               for (i in 1:nrow(referenceArea)){
                 referenceMatrix[referenceMatrix == referenceArea[i,1]] = 
                   1/referenceArea$Freq[i]                 
               }
             }
             )

    }
  }
  if (type == "Kneighbours"){
    # find neighbourhood
      Kneighbourhood = get.knnx(data = coordinates(simulations@locations),
               query = coordinates(simulations@locations)[referenceSD,],
               k = k)$nn.index      
    }
  
  # compute similarity
  similarity = numeric(length(SD))
  for (i in seq(along = SD)){
    switch(type,
             "equal" = {
               similarity[i] = length(intersect(SD[[i]], referenceSD))
               },
             "Kneighbours" = {
               similarity[i] = length(intersect(SD[[i]], Kneighbourhood))
             },
           "EarthMoversDistance" = {
             SDMatrix = simMatrix
             switch(class(simulations@locations),
                    "SpatialPixelsDataFrame" = {
                      SDMatrix = simMatrix
                      SDMatrix[!is.element(SDMatrix, SD[[i]]) & !is.na(SDMatrix)] = 0
                      SDMatrix[is.element(SDMatrix, SD[[i]])] = 1
                    },
                    "SpatialPolygridDataFrame" = {
                      SDMatrix = simMatrix
                      SDArea = as.data.frame(table(simMatrix[is.element(simMatrix, SD[[i]])]))
                      SDMatrix[!is.element(SDMatrix, SD[[i]]) & !is.na(SDMatrix)] = 0
                      for (j in 1:nrow(SDArea)){
                        SDMatrix[SDMatrix == SDArea[j,1]] = 
                          1/SDArea$Freq[j]                 
                      }
                    }
                    )
             similarity[i] = emd2d(referenceMatrix, SDMatrix,
                                   xdist = simulations@locations@grid@cellsize[1],
                                   ydist = simulations@locations@grid@cellsize[2], ...)
           }
    )    
  }
  return(similarity)
}  


# goodness on plumes for calibration and test

