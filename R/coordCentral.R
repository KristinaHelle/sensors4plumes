#############################################################
#                   coordCentral                            #
#############################################################
# kind of 'centroids' for SpatialLinesDataFrame (mainly for plotting and identification in optimiseSD_manual)
## coordinates of "corner" which is closest to mean coordinate of all corners (of each feature) 

coordCentral = function(x){# "corner" coordinate of line which is closest to centroid
  switch(class(x),
         "SpatialLinesDataFrame" = {
           coordAll = coordinates(x)
           coordCentral = matrix(nrow = length(x), ncol = 2)
           for (i in seq(along = x)){
             nCorners_i = sapply(X = coordAll[[i]], FUN = nrow)
             coordAll_i = matrix(nrow = sum(nCorners_i), ncol = 2)
             h = 0
             for (j in seq(along = coordAll[[i]])){
               coordAll_i[h + 1:nCorners_i[j],] = coordAll[[i]][[j]]  
               h = h + nCorners_i[j]
             }
             coordMean_i = matrix(apply(X = coordAll_i, FUN = mean, MARGIN = 2), ncol = 2)
             nearest_i = get.knnx(data = coordAll_i, query = coordMean_i, k = 1)
             coordCentral[i,] = coordAll_i[nearest_i$nn.index,]
           }
           coordCentral = as.data.frame(coordCentral)
           coordinates(coordCentral) = 1:2
           proj4string(coordCentral) = CRS(proj4string(x))
         },
         "SpatialIndexDataFrame" = {
           coordCentral = data.frame(y = length(x@index):1, x = 1)
           coordCentral = coordCentral[!duplicated(x@index),]
           coordCentral = coordCentral[unique(x@index),]
           coordinates(coordCentral) = ~ x + y
         },{
           coordCentral = as.data.frame(coordinates(x))
           coordinates(coordCentral) = 1:2
           proj4string(coordCentral) = CRS(proj4string(x))
         }
         )
  return(coordCentral)
}
