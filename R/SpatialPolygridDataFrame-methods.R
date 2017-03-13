###############################################################
#              SpatialPolygridDataFrame-methods               #
###############################################################
# coerce: SpatialGridDataFrame ->
#         -> SpatialGridDataFrame
#         SpatialPointsDataFrame ->
#         -> SpatialPointsDataFrame

# coordinates
# proj4string
# spplot
# cbind
# subsetSDF => subsetSDF
# areaSDF => areaSDF


setAs("SpatialGridDataFrame", "SpatialPolygridDataFrame",
      function(from){
        to = new("SpatialPolygridDataFrame", 
                 data = from@data, 
                 index = seq(length = nrow(from@data)),
                 grid = from@grid,
                 bbox = bbox(from),
                 proj4string = from@proj4string)   
        return(to)
      })


setAs("SpatialPolygridDataFrame", "SpatialGridDataFrame",
      function(from){
        to = polygrid2grid(
          obj = from, 
          zcol = names(from@data),  
          returnSGDF = TRUE
        )
        return(to)
      })


setAs("SpatialPointsDataFrame", "SpatialPolygridDataFrame",
      function(from){
        to = points2polygrid(
          points = SpatialPoints(from@coords, proj4string = from@proj4string),
          data = from@data)   
        return(to)
      })

setAs("SpatialPolygridDataFrame", "SpatialPointsDataFrame",
      function(from){
        to = SpatialPointsDataFrame(
          coords = coordinates(from),
          data = from@data,
          proj4string = CRS(proj4string(from))
        )
        return(to)
      })

setMethod("coordinates", signature(obj = "SpatialPolygridDataFrame"),
          function(obj){ 
            cells = SpatialGrid(obj@grid, obj@proj4string)
            fullgrid(cells) = FALSE
            cells = as(cells, "data.frame")            
            #coordinates_grid = coordinates(obj@grid)
            coord1 = split(cells[,1], obj@index)
            coord2 = split(cells[,2], obj@index)         
            coordinates = matrix(
              c(sapply(coord1, mean, na.rm = TRUE),
                sapply(coord2, mean, na.rm = TRUE)), 
              ncol = 2, byrow = FALSE)
            colnames(coordinates) = names(cells)
            return(coordinates)
          }
)

setMethod("proj4string", signature(obj = "SpatialPolygridDataFrame"),
          function(obj){ 
            proj4string = obj@proj4string@projargs
            return(proj4string)
          }
)

setMethod("bbox", signature(obj = "SpatialPolygridDataFrame"),
          function(obj){ 
            bbox = obj@bbox
            return(bbox)
          }
)

setMethod("spplot", signature(obj = "SpatialPolygridDataFrame"),
          function(obj, geoTiffPath, zcol = names(obj@data), plot = TRUE, returnSGDF = FALSE, ...){ 
            if (sum(obj@grid@cells.dim) > 2){
              result = polygrid2grid(
                obj = obj, 
                geoTiffPath = geoTiffPath, 
                zcol = zcol,  
                returnSGDF = TRUE)
              print(spplot(result, zcol = zcol,...))       
              if(!missing(geoTiffPath)){
                writeGDAL(dataset = result, drivername = "GTiff", fname = paste(geoTiffPath, ".tif", sep = ""))
              }
            } else {
              stop("Cannot plot single grid cell, no grid returned, no geoTiff generated!")
            }
            if(returnSGDF){
              out = list()
              out[["spplot"]] = spplot(result, zcol = zcol, ...)
              out[["grid"]] = result
            } else {
              out = spplot(result, zcol = zcol,...)
            }
            return(out)
          }
)

setMethod("length", signature(x = "SpatialPolygridDataFrame"),
          function(x){ 
            length = nrow(x@data)
          }  
)

cbind.SpatialPolygridDataFrame = function(...){ 
  dots = list(...)
  n = length(dots)
  if (n >= 1){
    newObj = dots[[1]]   
  }else{
    stop("No object(s) given! Returns NULL.")
    newObj = NULL
  }
  if (n >= 2){
    for (l in 2:n){
      if(!identical(newObj@grid, dots[[l]]@grid)){
        stop("SpatialPolygridDataFrames cannot be combined as their 'grid' differ.")
      }else{
        if(!identical(newObj@proj4string, dots[[l]]@proj4string)){
          stop("SpatialPolygridDataFrames cannot be combined as their 'proj4string' differ.")
        }else{  
          if(!identical(newObj@index, dots[[l]]@index)){
            stop("SpatialPolygridDataFrames cannot be combined as their 'index' differ.")
          }else{
            newObj = new("SpatialPolygridDataFrame",
                         grid = newObj@grid,
                         index = newObj@index,
                         proj4string = newObj@proj4string,
                         bbox = newObj@bbox,
                         data = cbind(newObj@data, dots[[l]]@data)
            )
          }                
        }
      }   
    }
  }
  return(newObj)
}
