################################################################
#       loadSimulations                                        #
################################################################

loadSimulations = function(
  basicPath = ".", # path to basic directory with plume simulations files OR SpatialDataFrame with all plumes and kinds (name convention of data columns as else for file names)
  readBy = NA, # "scan" or "raster"  
  multilayer = "kinds", # if there are multiple layers/columns in the files, these mean different kind (else use "plumes")
  # the following parameters can be used if only part of the locations should be used, if the parameters are missing, the full region is used
  region, # logical, representing all locations in original data (if bBox used, it must refer to locations in bBox only); keep locations where TRUE; can be combined with bBox
  bBox, # xmin, xmax, ymin, ymax; bounding box if only part of the data should be used (works only for raster data)
  nameSave = NA, # path to directory where result shall be saved (raster needs a file to be saved, tmpfile does not work)
  overwrite = FALSE, # may output files be overwritten?
  ... # parametes to be forwarded to read.table and scan
){
  # get filenames and check if they fulfil the requirements
  paths0 = list.files(basicPath, recursive = FALSE)
  if (length(paths0) == 0){
    stop("No files in 'basicPath'.")
  }
  paths1 = sort(paths0)
  # split by suffix
  suffix0 = strsplit(paths1, "\\.")
  suffix0i = sapply(FUN = length, suffix0) # must all be 2
  if (any(suffix0i != 2)){
    warning("Not all files in the basicPath folder have a correct name, either their name contains '.' (in addition to the '.' before the suffix) or they have no suffix; these files are ignored.")
    paths1a = paths1[suffix0i == 2]
    paths1 = paths1a
    suffix0 = strsplit(paths1, "\\.")
  } 
  suffix1 = paste0(".", mapply(FUN = "[[", suffix0, 
                               sapply(FUN = length, X = suffix0)))
  paths2 = unlist(mapply(strsplit, paths1, suffix1))
  paths3 = strsplit(paths2, "_")  
  paths3i = sapply(FUN = length, paths3)
  kindsFile = FALSE
  if (all(paths3i == 2)){
    kindsFile = TRUE
    paths4 = cbind(paths1,
                   matrix(unlist(paths3), byrow = TRUE, ncol = 2))
  } else {
    if (all(paths3i == 1)){
      kindsFile = FALSE
      paths4 = cbind(paths1, unlist(paths3))
    } else {
      warning("Some file names do not fit the required pattern. Files are regardes as each containing a different plume and all the same kind(s) of plumes.")
      kindsFile = FALSE
      paths4 = cbind(paths1, paths2)
    }
  }
  
  # test suffix
  suffix2 = sort(unique(suffix1))
  knownSuffix = logical(4)
  names(knownSuffix) = c(".txt", ".csv", ".tif", ".grd")
  if (identical(suffix2, c(".grd", ".gri"))){
    if (all(paths2[suffix1  == ".grd"] == paths2[suffix1  == ".gri"])){
      paths4 = paths4[suffix1 == ".grd", , drop = FALSE]    
      knownSuffix[4] = TRUE
    } else {
      stop("There must be a '.grd' and a '.gri' file for each name.")
    }
  } else {
    if (length(suffix2) > 1){
      stop("Files have different suffixes, they must all have the same (or be pairs of '.grd' and '.gri').")
    } else {
      knownSuffix[1:3] = is.element(c(".txt", ".csv", ".tif"), suffix2)      
    }   
  }
  # method to read files, either taken from suffixes or parameter
  if (!any(knownSuffix) & missing(readBy)){
    stop("The suffix of the files is none of '.txt', '.csv', '.tif', '.grd' therefore the method to read the files must be provided as 'readBy'.")
  } else {
    if (!missing(readBy)){
            switch (readBy,
              "scan" = {readByScan = TRUE},
              "raster" = {readByScan = FALSE},
              stop("'readBy' must be 'scan' or 'raster'."))    
    } else {
      if (any(knownSuffix[1:2])){
        readByScan = TRUE
      } else {
        readByScan = FALSE
      }
    }
  }
  
  # test if names (split by kind) are unique and the same names occur for each kind 
  # determine number of kinds and of plumes
  if (kindsFile){
    paths5 = split(paths4[,2], paths4[,3])
    for (i in seq(along = paths5)){
      if (length(unique(paths5[[i]])) < length(paths5[[i]])){
        stop(paste0("Names of plumes of kind ", names(paths5)[i], " are not unique."))
      }
      if (!identical(paths5[[i]], paths5[[1]])){
        stop(paste0("The plume names for kind ", names(paths5)[i], " differ from those for kind ", names(paths5)[1], "."))
      }
    }
    nP = length(paths5[[1]])
    nK = length(paths5)
  } else {
    if (length(unique(paths4[,2])) < length(paths4[,2])){
      stop("Names are not unique.")
    }
    nP = nrow(paths4)
    nK = 1
  }



  # load example file and check for multicolumns, determine size
  if (!readByScan){
    data1 = brick(paste0(basicPath, "/", paths4[1,1]), ...)
    nLayers = nlayers(data1)    
    if (!missing(bBox)){
      data1 = crop(data1, bBox) 
    }
    nL = ncell(data1)
  } else {
    data1 = read.table(file = paste0(basicPath, "/", paths4[1,1]), ...)
    nLayers = ncol(data1)   
    nL = nrow(data1)
  }
  if (!missing(region)){
    if (length(region) != nL | !is.logical(region)){
      warning("'region' must be logical and it's length must fit number of locations (after cropping to 'bBox'); this is not the case, therefore 'region' is ignored, all original locations are used.")
      region = NULL
    } else {
      print("Extract data for the chosen region.")
      nL = sum(region == TRUE) 
    }
  } else {
    region = NULL
  }
  
  if (nLayers > 1){
    if (ncol(paths4) > 2){
      warning("File names indicate that files contain different kinds of values for the plumes, therefore only first layer/column of each file is used.")
      multilayer == NA  
    }
    if (multilayer == "kinds"){
      nK = nLayers
    }
    if (multilayer == "plumes"){
      nP = nP * nLayers
    }
    if (!is.element(multilayer, c("kinds", "plumes", NA))){
      warning("Invalid value for 'multilayer'; only first layer/column of each file is used.")
      multilayer = NA
    }
  } else {
    multilayer = NA
  }    
  
  print(paste("The data consists of ", nK, " different type(s) of values of ", nP, " plumes in ", nL, " locations.", sep = ""))
  
  
  if (readByScan){
    out = loadSimulations_scan(
      basicPath = basicPath,
      filePaths = paths4,
      region = region,
      nameSave = nameSave,
      overwrite = overwrite,
      nP = nP,
      nL = nL,
      nK = nK,
      multilayer = multilayer, 
      ...)    
  } else {
    out = loadSimulations_raster(
      basicPath = basicPath,
      filePaths = paths4,
      region = region,
      bBox = bBox,
      nameSave = nameSave,
      overwrite = overwrite,
      nP = nP,
      nL = nL,
      nK = nK,
      multilayer = multilayer,
      ...)
  }
  return(out)
}  

#setMethod("loadSimulations", signature(basicPath = "character"), loadSimulations)
