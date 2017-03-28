#########################################################################
# simulationsApply                                                      #
#########################################################################

# add warning if nout wrong in all cases (if it does not matter, use warning; else use stop)

simulationsApply = function(
  simulations, # "Simulations" or "Raster*"
  locations = 1:nLocations(simulations), # "integer", indices of simulations@locations to be used
  plumes = 1:nPlumes(simulations), # "integer", indices of simulations@plumes to be used
  kinds = 1:nKinds(simulations), # 'integer' or 'character'; layers to be kept
  fun = NA, # no MARGIN, to be applied to full set of selected data
  fun_p = NA,# MARGIN = p
  fun_l = NA,# MARGIN = l
  fun_pl = NA, # MARGIN = c(p,l); function(x, weight = "source", nout){weight * x[,1]}
  fun_Rp = NA, # to be applied to result of fun_p
  fun_Rl = NA, # to be applied to result of fun_l
  fun_Rpl = NA,# to be applied to result of fun_pl
  fun_Rpl_cellStats = NA, # as fun_Rpl, to be used if not all data in memory; c("sum", "mean", "min", "max", "sd", "skew", "rms")
  nameSave = "simulationsApply", # filename, if results don't fit into memory: simulationsApply_global.grd etc.; if FALSE: nothing saved to file, stop if it would be necessary
  overwrite = FALSE,
  chunksize = 1e+7, #raster:::.chunksize(),
  keepSubset = FALSE,
  ...
){
  #----------------- adapt to class of simulations ----------------------
  if (is.element(class(simulations), c("RasterLayer", "RasterStack", "RasterBrick"))){
    data = simulations
  } else {
    if (is(simulations, "Simulations")){
      data = simulations@values
    } else {
      stop("'simulations' must be of class 'Simulations' or of a 'Raster*' class.")
    }
  }

  #--------------- determine subset size -----------------------------
  # if not all values to be used, delete superflous
  data = subset(data, kinds)
  nLay = nlayers(data)

  # determine plumes to be used
  nP = ncol(data)
  if (!identical(plumes, 1:nP)){
    plumesIn = plumes > 0 & plumes <= nP
    plumes[!plumesIn] = NA
    if (any(is.na(plumes))){
      stop("'plumes' out of bounds or contains 'NA'.")
    }
    nPl = length(plumes)
    isPlumes = TRUE
  } else {
    nPl = nP
    nPlu = nP
    plumes = 1:nPl
    plumesU = 1:nPlu
    plumesIndex = 1:nP
    isPlumes = FALSE
  }

  # determine locations to be used
  nL = nrow(data)
  if (!identical(locations, 1:nL)){
    locationsIn = locations > 0 & locations <= nL
    locations[!locationsIn] = NA
    if (any(is.na(locations))){
      stop("'locations' out of bounds or contains 'NA'.")
    }
    nLo = length(locations) # final length, including repetitions
    isLocations = TRUE
  }else{
    nLo = nL
    nLoc = nL
    locations = 1:nLo
    locationsU = 1:nLoc
    locationsIndex = 1:nL
    isLocations = FALSE
  }

  #--------------- check functions to be applied -----------------------
  functions = list()
  functionValid = logical(8)
  names(functionValid) = c("fun", "fun_l", "fun_p", "fun_pl", "fun_Rl", "fun_Rp", "fun_Rpl", "fun_Rpl_cellStats")
  if (is.function(fun)){
    fun_ = replaceDefault(fun, type = "fun.simulationsApply")
    functionValid["fun"] = fun_[[2]]
    functions[["fun"]] = fun_[[1]]
  }
  if (is.function(fun_l)){
    fun_l_ = replaceDefault(fun_l, type = "fun.simulationsApply")
    functionValid["fun_l"] = fun_l_[[2]]
    functions[["fun_l"]] = fun_l_[[1]]
  }
  if (is.function(fun_p)){
    fun_p_ = replaceDefault(fun_p, type = "fun.simulationsApply")
    functionValid["fun_p"] = fun_p_[[2]]
    functions[["fun_p"]] = fun_p_[[1]]
  }
  if (is.function(fun_pl)){
    fun_pl_ = replaceDefault(fun_pl, type = "fun.simulationsApply")
    functionValid["fun_pl"] = fun_pl_[[2]]
    functions[["fun_pl"]] = fun_pl_[[1]]
  }
  if (is.function(fun_Rl)){
    if (is(simulations, "Simulations")){
      fun_Rl_ = replaceDefault(fun_Rl,
                               newDefaults = list(weight = simulations@locations@data[locations,, drop = FALSE]),
                               type = "funR.simulationsApply")
    } else {
      fun_Rl_ = replaceDefault(fun_Rl, #-> fun must have 'weight', default is not changed
                               type = "funR.simulationsApply")
    }
      functionValid["fun_Rl"] = fun_Rl_[[2]]
      functions[["fun_Rl"]] = fun_Rl_[[1]]
  }
  if (is.function(fun_Rp)){
    if (is(simulations, "Simulations")){
      fun_Rp_ = replaceDefault(fun_Rp,
                               newDefaults = list(weight = simulations@plumes[plumes,, drop = FALSE]),
                               type = "funR.simulationsApply")
    } else {
      fun_Rp_ = replaceDefault(fun_Rp, #-> fun must have default for weight!
                               type = "funR.simulationsApply")
    }
      functionValid[["fun_Rp"]] = fun_Rp_[[2]]
      functions[["fun_Rp"]] = fun_Rp_[[1]]
  }
  if (is.function(fun_Rpl)){
    if (is(simulations, "Simulations")){
      fun_Rpl_ = replaceDefault(fun_Rpl,
                                newDefaults = list(
                                  weight_l = simulations@locations@data[locations,, drop = FALSE],
                                  weight_p = simulations@plumes[plumes,, drop = FALSE]),
                                type = "funRR.simulationsApply")
    } else {
      fun_Rpl_ = replaceDefault(fun_Rpl, type = "funRR.simulationsApply")
    }
    functionValid["fun_Rpl"] = fun_Rpl_[[2]]
    functions[["fun_Rpl"]] = fun_Rpl_[[1]]
  }
  if (!is.na(fun_Rpl_cellStats)){
    if (is.element(fun_Rpl_cellStats,c("sum", "mean", "min", "max", "sd", "skew", "rms"))){
      functionValid["fun_Rpl_cellStats"] = TRUE
    } else {
      warning("'fun_Rpl_cellStats' cannot be used, it has to be one of the strings: c('sum', 'mean', 'min', 'max', 'sd', 'skew', 'rms').")
    }
  }

  if (!all(functionValid[names(functions)])){
    warning(paste0("Some of the defined functions are invalid (e.g. because of missing parameters or extra parameters without default: ",
                   names(functions)[!functionValid[names(functions)]], ". Their result cannot be computed."))
  }
  # delete multiple indices and sort
  ## (must be done AFTER checkFun(fun_Rl) as there the original locations[except NA] are needed)
  if (isLocations){
    locationsTable = table(locations)
    locationsRank = rank(locations)
    locationsU = sort(unique(locations))
    isLocUnique = identical(locations, locationsU)
    nLoc = length(locationsU) # length of unique
    if (nLoc > 0){
      locationsIndex = unlist(mapply(rep, 1:nLoc, locationsTable))[locationsRank]
    } else {
      locationsIndex = integer(0)
    }
  } else {
    isLocUnique = TRUE
  }

  if (isPlumes){# delete multiple indices and sort
    plumesTable = table(plumes)
    plumesRank = rank(plumes)
    plumesU = sort(unique(plumes))
    isPluUnique = identical(plumes, plumesU)
    nPlu = length(plumesU) # length of unique
    if (nPlu > 0){
      plumesIndex = unlist(mapply(rep, 1:nPlu, plumesTable))[plumesRank]
    } else {
      plumesIndex = integer(0)
    }
  } else {
    isPluUnique = TRUE
  }

  #------------------ test which data can be held in memory -----------------------
  ## all data
  bs = blockSize(data, minblocks = 1, chunksize = chunksize)

  ## all selected data
  bs_subset = blockSize(raster(nrow = nLo, ncol = nPl), n = nLay, minblocks = 1, chunksize = chunksize)
  ## all results of functions
  if (functionValid["fun"]){
    bs_fun = blockSize(raster(nrow = 1, ncol = 1), n = eval(formals(functions[["fun"]])[["nout"]]),
                       minblocks = 1, chunksize = chunksize)
    if (bs_fun$n > 1){
      warning("Result of 'fun' too big to keep in memory, not returned.")
      functionValid["fun"] = FALSE
    }
  }
  if (functionValid["fun_pl"]){
    bs_fun_pl = blockSize(raster(nrow = nLo, ncol = nPl),
                          n = eval(formals(functions[["fun_pl"]])[["nout"]]),
                          minblocks = 1, chunksize = chunksize)
    if (bs_fun_pl$n > 1){
      if (nameSave == FALSE){
        warning("Result of 'fun_pl' too big to keep in memory, not returned.")
        functionValid[["fun_pl"]] = FALSE
      }else{
        rasterName_fun_pl = paste(nameSave, "_locationsplumes.grd", sep = "")
        warning(paste0("Result of 'fun_pl' too big to keep in memory,
                      saved at '", rasterName_fun_pl, "'."))
      }
    }
  }
  if (functionValid["fun_p"]){
    bs_fun_p = blockSize(raster(nrow = nPl, ncol = 1),
                         n = eval(formals(functions[["fun_p"]])[["nout"]]),
                         minblocks = 1, chunksize = chunksize)
    if (bs_fun_p$n > 1){
      warning("Result of 'fun_p' too big to keep in memory, not returned.")
      functionValid["fun_p"] = FALSE
    }
  }
  if (functionValid["fun_l"]){
    bs_fun_l = blockSize(raster(nrow = nLo, ncol = 1),
                         n = eval(formals(functions[["fun_l"]])[["nout"]]),
                         minblocks = 1, chunksize = chunksize)
    if (bs_fun_l$n > 1){
      warning("Result of 'fun_l' too big to keep in memory, not returned.")
      functionValid["fun_l"] = FALSE
    }
  }
  if (functionValid["fun_Rp"]){
    if (functionValid["fun_p"]){
      bs_fun_Rp = blockSize(raster(nrow = 1, ncol = 1),
                            n = eval(formals(functions[["fun_Rp"]])[["nout"]]),
                            minblocks = 1, chunksize = chunksize)
      if (bs_fun_Rp$n > 1){
        warning("Result of 'fun_Rp' too big to keep in memory, not returned.")
        functionValid["fun_Rp"] = FALSE
      }
    } else {
      warning("'fun_Rp' is to be applied to the results of 'fun_p', as 'fun_p' is missing or cannot be applied, no results of 'fun_Rp' returned.'")
      functionValid["fun_Rp"] = FALSE
    }
  }
  if (functionValid["fun_Rl"]){
    if (functionValid["fun_l"]){
      bs_fun_Rl = blockSize(raster(nrow = 1, ncol = 1),
                            n = eval(formals(functions[["fun_Rl"]])[["nout"]]),
                            minblocks = 1, chunksize = chunksize)
      if (bs_fun_Rl$n > 1){
        warning("Result of 'fun_Rl' too big to keep in memory, not returned.")
        functionValid["fun_Rl"] = FALSE
      }
    } else {
      warning("'fun_Rl' is to be applied to the results of 'fun_l', as 'fun_l' is missing or cannot be applied, no results of 'fun_Rl' returned.'")
      functionValid["fun_Rl"] = FALSE
    }
  }
  if (functionValid["fun_Rpl"]){
    if (functionValid["fun_pl"]){
      bs_fun_Rpl = blockSize(raster(nrow = 1, ncol = 1),
                             n = eval(formals(functions[["fun_Rpl"]])[["nout"]]),
                             minblocks = 1, chunksize = chunksize)
      if (bs_fun_pl$n > 1){
        warning("Result of 'fun_pl' not in memory, therefore 'fun_Rpl' cannot be applied.")
        functionValid["fun_Rpl"] = FALSE
      }
      if (bs_fun_Rpl$n > 1){
        warning("Result of 'fun_Rpl' too big to keep in memory, not returned.")
        functionValid["fun_Rpl"] = FALSE
      }
    } else {
      warning("'fun_Rpl' is to be applied to the results of 'fun_pl', as 'fun_pl' is missing or cannot be applied, no results of 'fun_Rpl' returned.'")
      functionValid["fun_Rpl"] = FALSE
    }
  }

  #--------------------- extract data (and apply functions) ---------------------------------
  result = list()

  if (bs_subset$n == 1){ # all extracted data can be loaded into memory
    result_subset = array(dim = c(nLoc, nPlu, nLay))
    h = 0
    for (i in 1:bs$n){
      # rows of this block
      locations_i = bs$row[i] - 1 + 1:bs$nrows[i]
      if (isLocations){
        which_locations_i = is.element(locationsU, locations_i)
        locations_i = locationsU[which_locations_i]
      }
      locations_i = locations_i - (bs$row[i] - 1)
      nLoc_i = length(locations_i)
      if (nLoc_i > 0){
        data_i = getValues(data, row = bs$row[i], nrows = bs$nrows[i])
        in_i_array = in_i_array = aperm(array(data_i, dim = c(nP, bs$nrows[i], nLay)), c(2,1,3)) # plumes, locations, values
        result_subset[h + 1:nLoc_i,,] = in_i_array[locations_i, plumesU,, drop = FALSE]
      }
      h = h + nLoc_i
    }
    # restore original order and multiplicity of locations and plumes
    if(isLocations){
      result_subset = result_subset[locationsIndex,,,drop = FALSE]
    }
    if(isPlumes){
       result_subset = result_subset[,plumesIndex,,drop = FALSE]
    }

    # apply functions
    if (functionValid["fun"]){
      result[["result_global"]] = functions[["fun"]](x = result_subset)
    }

    if (functionValid["fun_l"]){
      result[["result_locations"]] = apply(X = result_subset, FUN = functions[["fun_l"]], MARGIN = 1)
      if (is.null(dim(result[["result_locations"]]))){
        result[["result_locations"]] = as.matrix(result[["result_locations"]])
      } else {
        result[["result_locations"]] = t(result[["result_locations"]])
      }

      if (functionValid["fun_Rl"]){
        result[["result_global_locations"]] = functions[["fun_Rl"]](x = result[["result_locations"]])
      }
    }

    if (functionValid["fun_p"]){
      result[["result_plumes"]] = apply(X = result_subset, FUN = functions[["fun_p"]], MARGIN = 2)
      if (is.null(dim(result[["result_plumes"]]))){
        result[["result_plumes"]] = as.matrix(result[["result_plumes"]])
      } else {
        result[["result_plumes"]] = t(result[["result_plumes"]])
      }

      if (functionValid["fun_Rp"]){
        result[["result_global_plumes"]] = functions[["fun_Rp"]](x = result[["result_plumes"]])
      }
    }

    if (functionValid["fun_pl"]){
      if(bs_fun_pl$n <= 1){
        result[["result_locationsplumes"]] = apply(X = result_subset, FUN = functions[["fun_pl"]],
                                                   MARGIN = c(1,2))
        if (class(result[["result_locationsplumes"]]) != "RasterBrick"){
          if (length(dim(result[["result_locationsplumes"]])) == 2){
            dim(result[["result_locationsplumes"]]) = c(dim(result[["result_locationsplumes"]]), 1)
          } else {
            result[["result_locationsplumes"]] = aperm(result[["result_locationsplumes"]], perm = c(2,3,1))
          }
          if (prod(dim(result[["result_locationsplumes"]])) > 0){
            result[["result_locationsplumes"]] = brick(result[["result_locationsplumes"]],
                                                       xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                                                       crs = "+init=epsg:4326")
          } else {
            warning("As 'result[['result_locationsplumes']]' has size 0, it is not transformed into a brick and 'fun_Rpl' cannot be applied.")
            functionValid["fun_Rpl"] = FALSE
          }
        }

        if (functionValid["fun_Rpl"]){
          result[["result_global_locationsplumes"]] =
            functions[["fun_Rpl"]](x = getValues(result[["result_locationsplumes"]]))
        }
      } else {
        if(isLocUnique){
          result[["result_locationsplumes"]] = brick(nrows = nLoc,
                                                     ncols = nPl,
                                                     nl = eval(formals(functions[["fun_pl"]])$nout), # this can be assumed to be >1, else it could be kept in memory
                                                     xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                                                     crs = "+init=epsg:4326")
          print(paste0("write raster to ", rasterName_fun_pl))
          result[["result_locationsplumes"]] = writeStart(result[["result_locationsplumes"]],
                                                          filename = rasterName_fun_pl,
                                                          overwrite = overwrite)
          for (i in 1:bs_fun_pl$n) {
            x_i = result_subset[bs_fun_pl$row[i] + 1:bs_fun_pl$nrows[i] - 1,,]     # load some rows
            x_i = aperm(x_i, perm = c(2,1,3))
            dim(x_i) = c(bs_fun_pl$nrows[i] * nPl, nLay)                             # turn into cells x layers
            out_i = t(apply(X = x_i, FUN = functions[["fun_pl"]], MARGIN = 1))                                                      # apply fun cell-wise
            true_nout = dim(out_i)[2]
            if (!(eval(formals(functions[["fun_pl"]])$nout) == true_nout)){
              warning(paste0("The argument 'nout' of 'fun_pl' must fit the actual length of output of 'fun_pl'
                         if applied to simulations@values[i,j,] for any i, i. For the chosen function nout = ",
                             true_nout, ". No result of 'fun_pl' and 'fun_Rpl' returned."))
              functionValid[c("fun_pl", "fun_Rpl")] = FALSE
              result[["result_locationsplumes"]] = writeStop(result[["result_locationsplumes"]])
              result[["result_locationsplumes"]] = NULL
              break
            }
            if(!isPluUnique){
              # repeat/reorder
              out_i = out_i[rep(plumesIndex, bs_fun_pl$nrows[i]),]
            }
            writeValues(result[["result_locationsplumes"]], out_i, bs_fun_pl$row[i])          # save result to output raster object
          }
          if (functionValid["fun_pl"]){
            result[["result_locationsplumes"]] = writeStop(result[["result_locationsplumes"]])
          }
        }else{# if !isLocUnique: create intermediate result first, later copy with  reorder/repeat
         result_locationsplumes = brick(nrows = nLoc,
                                        ncols = nPl,
                                        nl = eval(formals(fun_pl)$nout), # this can be assumed to be >1, else it could be kept in memory
                                        xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                                        crs = "+init=epsg:4326")
          warning("Intermediate results are saved in 'intermediateResult_locationsplumes.grd', the file is deleted in the end.")
          print("write raster to 'intermediateResult_locationsplumes.grd'")
          result_locationsplumes = writeStart(result_locationsplumes,
                                              filename = "intermediateResult_locationsplumes.grd",
                                              overwrite = overwrite)
          for (i in 1:bs_fun_pl$n) {
            x_i = result_subset[bs_fun_pl$row[i] + 1:bs_fun_pl$nrows[i] - 1,,]  # load some rows
            x_i = aperm(x_i, perm = c(2,1,3))
            dim(x_i) = c(bs_fun_pl$nrows[i] * nPl, nLay)                             # turn into cells x layers
            out_i = t(apply(X = x_i, FUN = fun_pl, MARGIN = 1))                                                      # apply fun cell-wise
            true_nout = dim(out_i)[2]
            if (!(eval(formals(fun_pl)$nout == true_nout))){
              warning(paste0("The argument 'nout' of 'fun_pl' must fit the actual length of output of 'fun_pl'
                         if applied to simulations@values[i,j,] for any i, i. For the chosen function nout = ",
                             true_nout, ". No result of 'fun_pl' and 'fun_Rpl' returned."))
              functionValid[c("fun_pl", "fun_Rpl")] = FALSE
              result_locationsplumes = writeStop(result_locationsplumes)
              result[["result_locationsplumes"]] = NULL
              break
            }
            if(!isPluUnique){
              # repeat/reorder
              out_i = out_i[rep(plumesIndex, bs_fun_pl$nrows[i]),]
            }
            writeValues(result_locationsplumes, out_i, bs_fun_pl$row[i])          # save result to output raster object
          }
          if (functionValid["fun_pl"]){
            result_locationsplumes = writeStop(result_locationsplumes)

            # create final result with correct order of locations
            result[["result_locationsplumes"]] = brick(nrows = nLoc,
                                                       ncols = nPlu,
                                                       nl = eval(formals(functions[["fun_pl"]])$nout), # this can be assumed to be >1, else it could be kept in memory
                                                       xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                                                       crs = "+init=epsg:4326")
            print(paste0("write raster to ", rasterName_fun_pl))
            result[["result_locationsplumes"]] = writeStart(result[["result_locationsplumes"]],
                                                            filename = rasterName_fun_pl,
                                                            overwrite = overwrite)
            for (i in seq(along = locationsIndex)){
              writeValues(result[["result_locationsplumes"]],
                          getValues(result_locationsplumes, row = locationsIndex[i]),
                          i)
            }
            result[["result_locationsplumes"]] = writeStop(result[["result_locationsplumes"]])
         }
         rm(result_locationsplumes)
         print("remove file 'intermediateResult_locationsplumes.grd'.")
         file.remove("intermediateResult_locationsplumes.grd")
         file.remove("intermediateResult_locationsplumes.gri")
        }
      }
    }
  }else{ # subset cannot be loaded at once
    if (functionValid["fun"]){
      warning("'fun' cannot be applied as this would require to load all selected data into memory at once.")
      functionValid["fun"] = FALSE
    }
    if (functionValid["fun_l"]){
      # define where to save the result
      result[["result_locations"]] = array(dim = c(nLoc, eval(formals(functions[["fun_l"]])$nout)))
      # load and process in chunks
      h = 0
      for (i in 1:bs$n){
        # rows of this block
        locations_i = bs$row[i] - 1 + 1:bs$nrows[i]
        if (isLocations){
          which_locations_i = is.element(locationsU, locations_i)
          locations_i = locationsU[which_locations_i]
        }
        locations_i = locations_i - (bs$row[i] - 1)
        nLoc_i = length(locations_i)
        if (nLoc_i > 0){
          data_i = getValues(data, row = bs$row[i], nrows = bs$nrows[i])
          in_i_array = aperm(array(data_i, dim = c(nP, bs$nrows[i], nLay)), c(2,1,3))
          result[["result_locations"]][h + 1:nLoc_i,] = t(apply(X = in_i_array[locations_i,plumes,,drop = FALSE],
                                                                FUN = functions[["fun_l"]], MARGIN = 1))
        }
        h = h + nLoc_i
      }
      if(isLocations){
        result[["result_locations"]] = result[["result_locations"]][locationsIndex,,drop=FALSE]
      }
      if (functionValid["fun_Rl"]){
        result[["result_global_locations"]] = functions[["fun_Rl"]](x = result[["result_locations"]])
      }
    }
    if (functionValid["fun_p"]){
      # define where to save the result
      result[["result_plumes"]] = array(dim = c(nPlu, eval(formals(functions[["fun_p"]])$nout)))
      # load and process in chunks
      bs_columns = blockSize(raster(nrow = nPlu, ncol = nLo), n = nLay, minblocks = 1, chunksize = chunksize)
      for (g in 1:bs_columns$n){
        # columns of this block
        plumes_g = bs_columns$row[g] + 1:bs_columns$nrows[g] - 1
        nPlu_g = length(plumes_g)

        # define where to save the current block of columns
        result_g = array(dim = c(nLoc, nPlu_g, nLay))
        # load and process in chunks
        h = 0
        for (i in 1:bs$n){
          # rows of this block
          locations_i = bs$row[i] - 1 + 1:bs$nrows[i]
          if (isLocations){
            which_locations_i = is.element(locationsU, locations_i)
            locations_i = locationsU[which_locations_i]
          }
          locations_i = locations_i - (bs$row[i] - 1)
          nLoc_i = length(locations_i)
          if (nLoc_i > 0){
            data_i = getValues(data, row = bs$row[i], nrows = bs$nrows[i])
            in_i_array = aperm(array(data_i, dim = c(nP, bs$nrows[i], nLay)), c(2,1,3))
            result_g[h + 1:nLoc_i,,] = in_i_array[locations_i,plumesU[plumes_g],,drop=FALSE]
          }
          h = h + nLoc_i
        }
        if(isLocations){
          result_g = result_g[locationsIndex,,,drop=FALSE]
        }
        result_plumes = apply(X =  result_g, FUN = functions[["fun_p"]], MARGIN = 2)
        if (is.vector(result_plumes)){
          if (!eval(formals(functions[["fun_p"]])[["nout"]]) == 1){
            warning("The argument 'nout' of 'fun_p' must fit the actual length of output of 'fun_p'
                    if applied to simulations@values[i,,] for any i. For the chosen function nout = 1.
                    No result of 'fun_p' returned.")
            functionValid["fun_p"] = FALSE
            result[["result_plumes"]] = NULL
            break
          } else {
            dim(result_plumes) = c(length(result_plumes), 1)
          }
        }else{
          true_nout = dim(result_plumes)[1]
          if (!(eval(formals(fun_p)[["nout"]]) == true_nout)){
            warning(paste("The argument 'nout' of 'fun_p' must fit the actual length of output of 'fun_p'
                          if applied to simulations@values[i,,] for any i. For the chosen function nout = ", true_nout, ".", sep = ""))
            functionValid["fun_p"] = FALSE
            result[["result_plumes"]] = NULL
            break
          }
        }
        if (functionValid["fun_p"]){
          result[["result_plumes"]][bs_columns$row[g] + 1:bs_columns$nrows[g] - 1,] = t(result_plumes)
        }
      }
      if (functionValid["fun_p"]){
        if(isPlumes){
          result[["result_plumes"]] = result[["result_plumes"]][plumesIndex,,drop=FALSE]
        }
        if (functionValid["fun_Rp"]){
          result[["result_global_plumes"]] = functions[["fun_Rp"]](x = result[["result_plumes"]])
        }
      }
    }
    if (functionValid["fun_pl"]){
      # define where to save the result
      if (bs_fun_pl$n <= 1){
        result[["result_locationsplumes"]] = array(dim = c(eval(formals(functions[["fun_pl"]])$nout), nLoc, nPlu))
      }else{
        # create raster files to save (intermediate) result
        if(isLocUnique){
          result[["result_locationsplumes"]] = brick(nrows = nLoc,
                                                     ncols = nPlu,
                                                     nl = eval(formals(functions[["fun_pl"]])$nout), # this can be assumed to be >1, else it could be kept in memory
                                                     xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                                                     crs = "+init=epsg:4326")
          print(paste0("write raster to ", rasterName_fun_pl))
          result[["result_locationsplumes"]] = writeStart(result[["result_locationsplumes"]],
                                                          filename = rasterName_fun_pl,
                                                          overwrite = overwrite)
        } else { # if !locUnique
          result_locationsplumes = brick(nrows = nLoc,
                                         ncols = nPl,
                                         nl = eval(formals(functions[["fun_pl"]])$nout), # this can be assumed to be >1, else it could be kept in memory
                                         xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                                         crs = "+init=epsg:4326")

          warning("Intermediate results are saved in 'intermediateResult_locationsplumes.grd', the file is deleted in the end.")
          print("write raster to 'intermediateResult_locationsplumes.grd'")
          result_locationsplumes = writeStart(result_locationsplumes,
                                              filename = "intermediateResult_locationsplumes.grd",
                                              overwrite = overwrite)

        }
      }
      # load and process in chunks, choose smallest necessary chunk size
      if (bs_subset$n < bs_fun_pl$n){
        bs_Fun_pl = bs_fun_pl
      }else{
        bs_Fun_pl = bs_subset
      }
      h = 0
      for (i in 1:bs_Fun_pl$n){
        # rows of this block
        locations_i = bs_Fun_pl$row[i] - 1 + 1:bs_Fun_pl$nrows[i]
        if (isLocations){
          which_locations_i = is.element(locationsU, locations_i)
          locations_i = locationsU[which_locations_i]
        }
        locations_i = locations_i - (bs_Fun_pl$row[i] - 1)
        nLoc_i = length(locations_i)
        if (nLoc_i > 0){
          data_i = getValues(data, row = bs_Fun_pl$row[i], nrows = bs_Fun_pl$nrows[i])
          if(is.vector(data_i)){
            data_i = matrix(data_i, ncol = 1)
          }
          in_i_array = aperm(array(data_i, dim = c(nP, bs_Fun_pl$nrows[i], nLay)), c(2,1,3)) # turn into loc x plu x val
          in_i_array_subset = in_i_array[locations_i, plumesU,,drop=FALSE]                   # subset
          if (bs_fun_pl$n <= 1){
            result[["result_locationsplumes"]][,h + 1:nLoc_i,] = apply(X = in_i_array_subset,
                                                                       FUN = functions[["fun_pl"]], MARGIN = c(1,2))
          }else{
            out_ii = apply(X = in_i_array_subset, FUN = functions[["fun_pl"]], MARGIN = c(1,2))
            if (length(dim(out_ii)) == 2){
              out_iii = array(dim = c(1, dim(out_ii)))
              out_iii[1,,] = out_ii
              out_ii = out_iii
            }
            out_i = aperm(out_ii, c(3,2,1))
            true_nout = dim(out_i)[3]
            if (eval(formals(fun_pl)$nout) != true_nout){
              warning(paste("The argument 'nout' of 'fun_pl' must fit the actual length of output of 'fun_pl'
                            if applied to simulations@values[i,j,] for any i, i. For the chosen function nout = ", true_nout, ".", sep = ""))
              functionValid["fun_pl"] = FALSE
              result[["result_locationsplumes"]] = writeStop(result[["result_locationsplumes"]])
              result[["result_locationsplumes"]] = NULL
              print(paste0("remove raster file ", rasterName_fun_pl))
              file.remove(rasterName_fun_pl)
              file.remove(paste0(strsplit(rasterName_fun_pl, ".grd")[[1]], ".gri"))
              break
            }
            if (functionValid["fun_pl"]){
              if (!isPluUnique){
                out_i = out_i[plumesIndex,,,drop = FALSE]
              }
              dim(out_i) = c(nLoc_i * nPl, eval(formals(fun_pl)$nout))              # turn into cells x layers
              if (isLocUnique){
                writeValues(result[["result_locationsplumes"]], out_i, h + 1)    # save result to output raster object
              } else {
                writeValues(result_locationsplumes, out_i, h +1)          # save result to output raster object
              }
            }
          }
        }
        h = h + nLoc_i
      }
      if (functionValid["fun_pl"]){
        if (isPluUnique) {
          result[["result_locationsplumes"]] = writeStop(result[["result_locationsplumes"]])
        } else {
          result_locationsplumes = writeStop(result_locationsplumes)
          result[["result_locationsplumes"]] = brick(nrows = nLo,
                                                     ncols = nPl,
                                                     nl = eval(formals(fun_pl)$nout), # this can be assumed to be >1, else it could be kept in memory
                                                     xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                                                     crs = "+init=epsg:4326")
          print(paste0("write raster to ", rasterName_fun_pl))
          result[["result_locationsplumes"]] = writeStart(result[["result_locationsplumes"]],
                                                          filename = rasterName_fun_pl,
                                                          overwrite = overwrite)

          for (i in seq(along = locationsIndex)){
            writeValues(result[["result_locationsplumes"]],
                        getValues(result_locationsplumes, row = locationsIndex[i]),
                        i)
          }
          result[["result_locationsplumes"]] = writeStop(result[["result_locationsplumes"]])
          print("remove raster 'intermediateResult_locationsplumes.grd'")
          rm(result_locationsplumes)
          file.remove("intermediateResult_locationsplumes.grd")
          file.remove("intermediateResult_locationsplumes.gri")
        }
      }
      if (functionValid["fun_pl"]){
        if (bs_fun_pl$n <= 1){
          if (isLocations){
            result[["result_locationsplumes"]] = result[["result_locationsplumes"]][,locationsIndex,,drop=FALSE]
          }
          if (isPlumes){
            result[["result_locationsplumes"]] = result[["result_locationsplumes"]][,,plumesIndex,drop=FALSE]
          }
          if (is.function(fun_Rpl)){
            result[["result_global_locationsplumes"]] =
              functions[["fun_Rpl"]](x = getValues(result[["result_locationsplumes"]]))
          }
        }
      }
    }
  }
  if (keepSubset){
    if (bs_subset$n == 1){
      result[["subset"]] = result_subset
    } else {
      warning("Subset too big to keep, only results of functions returned.")
    }
  }
  if (functionValid["fun_Rpl_cellStats"]){
    if (is.element("result_locationsplumes", names(result))){
      result[["cellStats_global_locationsplumes"]] = cellStats(x = result[["result_locationsplumes"]],
                                               stat = fun_Rpl_cellStats, asSample = FALSE)
    }
  }
  return(result)
}




