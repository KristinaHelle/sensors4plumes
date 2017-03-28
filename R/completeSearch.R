################################################################################
# complete search as function                                                  #
################################################################################
#! changes since "completeSearch_function.r" (should all have no effect on the optimisation itself):
## subsetting of matrix with drop = FALSE replaces extra code for dim = 1
## J is now big enough to search all SDs with UP TO n sensors
## saving improved: SDs with less than desired sensors that detect everything or fulfill the aim are recorded
completeSearch = function(Detectable,
                          #Detectable = detectable,                              # matrix of locations(rows) and plumes(columns) of 0 - plume here not detectable and 1; rows not necessarily
                          n = 5,                                                # number of sensors
                          lowerLimit = 1,                                       # only search for sensors that detect at least so many plumes
                          increaseLimit = FALSE,                                # if TRUE, the lower limit is increased to the current best cost + 1 every time the current best cost >= current lowerLimit; if FALSE lowerLimit is constant
                          iterations = NA,                                      # maximal number of iterations, if no value given: number of all possibilities
                          nameSave = NA                                         # path to save results, if it is a directory, it has to end with /
  ){
  out = list()
  Time = date()

  # ------------------- allocate some parameters ------------------------------ #
  ## to know which sensors have been tried and which ones have to
  cCandidates = list()
  cColumns = list()                                                             # current columns of cDetectable: given current sensors 1:(i-1) cDetectable to find the i-th sensor has these rows and columns [[i]]
  cRows = list()                                                                # as cColumns
  # --------------------------- initialise ----------------------------------- #
  i = 1

  ## make Detectable unique
  rownames(Detectable) = 1:nrow(Detectable)
  colnames(Detectable) = 1:ncol(Detectable)
  detectableCol = apply(FUN = any, X = Detectable, MARGIN = 2)
  M = sum(detectableCol)
  if(any(!detectableCol)){
    print(paste("Warning: only", M, "column(s) contain(s) TRUE, those without are ignored."))
    Detectable = Detectable[,detectableCol, drop = FALSE]
  }

  if(M < lowerLimit){# aim at detection of more plumes than exist
    warning(paste("There are only", M, "columns. The given limit must not exceed this number."))
    i = 0
    break
  }

  if (M < n){# aim to use more than one sensor per plume
    warning(paste0(M, " column(s) contain(s) TRUE. Therefore at most ", M, " sensor(s) are/is sufficient to detect all plumes. The desired number of sensors is set to ", M, "."))
    n = min(n, M)
  }


  cColumns[[i]] = colnames(Detectable)
  cRows[[i]] = rownames(Detectable)[!duplicated(Detectable)]
  cDetectable = Detectable[cRows[[i]], cColumns[[i]], drop = FALSE]                           # remaining plumes at remaining location units (equal locations merged)

  if(dim(cDetectable)[1] < n){
    print(paste("There are only", dim(cDetectable)[1], "different locations, so you cannot search for", n, "sensors."))
    i = 0
    break
  }

  # ------------------- allocate parameters ---------------------------------- #
  ## to know which sensors have been tried and which ones have to
  cSD = rep(0,n)                                                                # current sensors, given as indices of cCandidates                                                               # as cColumns
  cSDPlumes = integer(n)                                                        # for current sensors: how many plumes are exclusively detected by this sensor (not by the earlier ones)
  cNPlumesSD = integer(n)

  ## for reporting
  finalSD = matrix(nrow = 0, ncol = n)                                          # matrix: each row is a SD that fullfills the lowerLimit (original rowname)
  allSD = matrix(nrow = 0, ncol = n)                                            # for report: paste current cSDs  below each other
  nRnC = list()                                                                 # for report: size of cDetectable during all steps

  ## iterations (tested SDs)
  if(is.na(iterations)){
    J = sum(choose(dim(Detectable)[1], 1:n))                                    # stop after at most J iterations (number of all possible SDs of at most n sensors on the locations of Detectable)
  }else{
    J = iterations
  }
  j = 1


  ## determine where to put first sensor
  cNPlumes = apply(FUN = sum, X = cDetectable, MARGIN = 1)
  orderCNPlumes = order(cNPlumes, decreasing = TRUE)

  cCandidates[[i]] = data.frame(
                     name = as.character(rownames(cDetectable)[orderCNPlumes]),
                     nPlumes = cNPlumes[orderCNPlumes])
  nPlumesMaxMatrix = matrix(NA, ncol = nrow(cCandidates[[i]]), nrow = n)
  for(k in 1:n){
    nPlumesMaxMatrix[k,] = c(cCandidates[[i]]$nPlumes[k:length(cNPlumes)],rep(0, k-1))
  }
  cCandidates[[i]]$nPlumesMax = apply(FUN = sum, X = nPlumesMaxMatrix, MARGIN = 2)
  cCandidates[[i]]$use = cCandidates[[i]]$nPlumesMax >= lowerLimit
  out[["cCandidates"]] = cCandidates[[1]]

#! save.image(file = paste(nameSave, "initialWorkspace.Rdata", sep = ""))
  for (k in 1:(n+1)){
    nRnC[[k]] = matrix(nrow = 0, ncol = 2)
  }
  nRnC[[1]] = rbind(nRnC[[1]], c(length(cRows[[1]]), length(cColumns[[1]])))

  # ------------------------ optimisation ------------------------------------ #
  while(i > 0){
    if(sum(cCandidates[[i]]$use, na.rm = TRUE) > cSD[i]){
      cSD[i] = cSD[i] + 1
      cNPlumesSD[i] = cCandidates[[i]]$nPlumes[cSD[i]]
      allSD = rbind(allSD, cSD)

      if(i < n){
        if(sum(cNPlumesSD) >= M){
          print(paste("Already", sum(cSD > 0), "row(s) contain(s) TRUE for all", M, "columns."))
          thisSD = integer(n)
          for(k in 1:i){
            thisSD[k] = as.character(cCandidates[[k]]$name[cSD[k]])
          }
          finalSD = rbind(finalSD, thisSD)
          lowerLimit = M
          if (increaseLimit){
            lowerLimit = M + 1
          }
          if(!is.na(nameSave)){
            save(finalSD, file = paste0(nameSave, "_finalSD.Rdata", sep = ""))
          }
          i = 0
          break
        }
        if(!increaseLimit){
          if(sum(cNPlumesSD) >= lowerLimit){
            print(paste("Already", sum(cSD > 0), "row(s) contain(s) TRUE for", lowerLimit, " (given limit) columns."))
            thisSD = integer(n)
            for(k in 1:i){
              thisSD[k] = as.character(cCandidates[[k]]$name[cSD[k]])
            }
            finalSD = rbind(finalSD, thisSD)
            if(!is.na(nameSave)){
              save(finalSD, file = paste0(nameSave, "_finalSD.Rdata", sep = ""))
            }
            break
          }
        }


        ## make Detectable unique
        cDetectable = Detectable[cRows[[i]], cColumns[[i]], drop = FALSE]
        hereDetected = cDetectable[as.character(cCandidates[[i]]$name[cSD[i]]),, drop = FALSE]
        cColumns[[i + 1]] = colnames(cDetectable)[!hereDetected]
        cRows_i = setdiff(cRows[[i]], as.character(cCandidates[[i]]$name[1:cSD[i]]))
        cDetectable = cDetectable[cRows_i, cColumns[[i + 1]], drop = FALSE]
        cRows[[i + 1]] = cRows_i[!duplicated(cDetectable)]

        cDetectable =  Detectable[cRows[[i + 1]], cColumns[[i + 1]], drop = FALSE]

        nRnC[[i + 1]] = rbind(nRnC[[i + 1]], c(length(cColumns[[i + 1]]), length(cRows[[i + 1]])))

        ## determine where to put next level sensor for fixed cSD[i]
        cNPlumes = apply(FUN = sum, X = cDetectable, MARGIN = 1)
        if(all(dim(cDetectable) == c(0,0))){
          i = 0
          break
        }

        orderCNPlumes = order(cNPlumes, decreasing = TRUE)

        cCandidates[[i + 1]] = data.frame(
                     name = as.character(cRows[[i + 1]][orderCNPlumes]),
                     nPlumes = cNPlumes[orderCNPlumes])
        nCandidates = nrow(cCandidates[[i + 1]])
        nPlumesMaxMatrix = matrix(NA, ncol = nCandidates + n, nrow = n - i)
        if((n - i) > 1){
          for(k in seq(along = nPlumesMaxMatrix[,1])){
            nPlumesMaxMatrix[k,] = c(c(cCandidates[[i + 1]]$nPlumes, rep(0, n))[k:(nCandidates + n)],rep(0, k-1))
          }
          cCandidates[[i + 1]]$nPlumesMax = apply(FUN = sum, X = nPlumesMaxMatrix[,1:nCandidates, drop = FALSE], MARGIN = 2)
        }else{
          cCandidates[[i + 1]]$nPlumesMax = cCandidates[[i + 1]]$nPlumes
        }
        cCandidates[[i + 1]]$use = cCandidates[[i + 1]]$nPlumesMax + sum(cNPlumesSD[1:i]) >= lowerLimit


        if(sum(cCandidates[[i + 1]]$use, na.rm = TRUE) > 0){
          i = i + 1
        }
      }else{ # i == n
        thisSD = integer(n)
        for(k in 1:n){
          thisSD[k] = as.character(cCandidates[[k]]$name[cSD[k]])
        }
        nPlumesThisSD = sum(apply(FUN = any, Detectable[thisSD,,drop = FALSE], MARGIN = 2))
        if(nPlumesThisSD >= lowerLimit){
          finalSD = rbind(finalSD, thisSD)
          if(!is.na(nameSave)){
            save(finalSD, file = paste0(nameSave, "_finalSD.Rdata", sep = ""))
          }
          Time = c(Time, Success = date())
          print(paste0("Success [", date(), "]: ", paste(thisSD, collapse = " ")))

          if(increaseLimit){
            lowerLimit = nPlumesThisSD + 1
            if(lowerLimit > ncol(Detectable)){
              print(paste("Found a set of", n, "row(s) that contains TRUE for all columns."))
              i = 0
              break
            }
          }else{
            lowerLimit = nPlumesThisSD
          }
            cCandidates[[1]]$use = cCandidates[[1]]$nPlumesMax >= lowerLimit
          if(n > 1){
            for(k in 2:n){
              cCandidates[[k]]$use = cCandidates[[k]]$nPlumesMax + sum(cNPlumesSD[1:(k-1)]) >= lowerLimit
            }
          }
        }
      }
    }else{ # sum(cCandidates[[i]]$use) <= cSD[i]
      cSD[i:n] = 0
      cNPlumesSD[i:n] = 0
      for(k in i:n){
        cCandidates[[k]] = NULL
        cRows[[k]] = NULL
        cColumns[[k]] = NULL
      }
      i = i - 1
    }
    # save result each time first sensor changes
    if(i == 1){
      if(!is.na(nameSave)){
        save(allSD, nRnC, file =  paste0(nameSave, "_SD_nRnC_", cSD[1], ".Rdata", sep = ""))
      }
      allSD = matrix(nrow = 0, ncol = n)

#!     save(Time, finalSD, cCandidates, lowerLimit, cColumns, cRows, file = paste(nameSave, "currentWorkspace_", cSD[1], ".Rdata", sep = ""))
      #!save(nRnC, file =  paste(nameSave, "nRnC_", cSD[1], ".Rdata", sep = ""))
      nRnC = list()
      for (k in 1:(n+1)){
        nRnC[[k]] = matrix(nrow = 0, ncol = 2)
      }
      nRnC[[1]] = rbind(nRnC[[1]], c(length(cRows[[1]]), length(cColumns[[1]])))
      Time = c(Time, Loop = date())
    }
    j = j + 1
    if (j > J){
      print(paste("Break because maximal number of sampling designs to be tested (",J,") is reached.", sep = ""))
      i = 0
      break
    }
  }
  # --------------------- output --------------------------------------------- #
  Time = c(Time, End = date())
#!  save.image(file = paste(nameSave, "finalWorkspace_", cSD[1], ".Rdata", sep = ""))
  if(!is.na(nameSave)){
    save.image(file = paste0(nameSave, "_finalWorkspace.Rdata", sep = ""))
  }
  out[["finalSD"]] = finalSD
  out[["lowerLimit"]] = lowerLimit
  out[["Time"]] = Time

  return(out)
}

