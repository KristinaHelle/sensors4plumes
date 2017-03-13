###########################################################
#        test optimiseSD_global                           #
###########################################################

# OK errors
# OK logical / 0-1 have same results 
# OK invalid values -> early stop
## OK choice of layer via 'detectable' (character or integer)
## OK numbers > 1
## OK NA
## OK NaN

# OK (final) detectable has 0 rows/columns -> stop
# OK (final) detectable has 1 column or row -> solve

# OK searching for all SDs
# OK find minimal number of sensors required to detect all plumes -> one optimal SD
# OK searching for one optimal SD (may have more sensors than required to detect all plumes)

# OK locAll/locFix -> 
## OK previous deletion of plumes/locations
## OK back-transformation (internal indices to total location indices) as expected
## OK reconstruction (multiple optima) 

# OK aimNumber -> size of SD (if not reduced)
# OK missing, locationsInitial given -> length of locInit defines aimNumber

# OK maxIterations < choose(locAll, aimNumber) -> stop early, return reasonable (not necessary optimal) result

# OK nameSave = NA -> no files generated
# OK nameSave != NA -> intermediate resuls saved there 
# (useful for reconstruction if break)


# verbatim TRUE -> relevant information given to follow progress
# verbatim FALSE -> warnings for unexpected results; errors

# OK completeSearch works as expected (paper test case)
# OK realistic example (radioactivePlumes_area)

# OK call via optimiseSD


#-------- data -------------------- #
data(SimulationsSmall)
SimulationsSmall@values$detect = getValues(SimulationsSmall@values)[,1] > getValues(SimulationsSmall@values)[,2]

data(radioactivePlumes_area)
radioactivePlumes_area@values$detectable = calc(
  radioactivePlumes_area@values$maxdose,
  fun = function(x){x >= 1e-7})
locAll_l = rep(seq(0, 49, 7), times = 8) + rep(seq(0, 49, 7), each = 8) * 50 + 1
# spplot(subsetSDF(radioactivePlumes_area@locations, locations = locAll_l))  

# test case from paper
detectionMatrix = matrix(c(0,1,0,0,0,1,
                          1,0,1,0,0,0,
                          0,1,1,1,0,0,
                          1,0,0,0,0,0,
                          0,0,0,1,1,0,
                          1,1,0,0,0,0,
                          0,0,1,0,0,1), 
                        nrow = 7, ncol = 6, byrow = TRUE)
## to test difference between logical and 0/1
detectionMatrixLogic = matrix(as.logical(detectionMatrix), nrow = 7, ncol = 6) 
detectionMatrixError = array(dim = c(7,6,3))
## to test errors in data
for (i in 1:3){
  detectionMatrixError[,,i] = detectionMatrix  
}
detectionMatrixError[3,2,1] = NA
detectionMatrixError[4,1,2] = 10
detectionMatrixError[6,5,3] = NaN

detectionRaster = raster(x = detectionMatrix,
             xmn = -90, xmx = 90, ymn = -90, ymx = 90,
             crs = "+init=epsg:4326")
detectionRasterLogic = raster(x = detectionMatrixLogic,
                              xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                              crs = "+init=epsg:4326") 
dataType(detectionRasterLogic) = "LOG1S"
detectionRasterErrors = brick(x = detectionMatrixError,
                         xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                         crs = "+init=epsg:4326")

completeExample = Simulations( 
  locations = SimulationsSmall@locations[1:7,],
  plumes = SimulationsSmall@plumes[1:6,],
  values = stack(detectionRaster, detectionRasterLogic, detectionRasterErrors))

# test case: similar to paper, with multiple locations
detectionMatrix2 = matrix(c(
  0,1,0,0,0,1,# 1
  0,1,0,0,0,1,# 1 
  1,0,1,0,0,0,# 2
  0,1,1,1,0,0,# 3
  0,0,0,1,1,0,# 5 
  1,0,0,0,0,0,# 4
  0,0,0,1,1,0,# 5
  1,1,0,0,0,0,# 6
  0,0,1,0,0,1 # 7
  ), nrow = 9, ncol = 6, byrow = TRUE)

completeExample2 = Simulations( 
  locations = SimulationsSmall@locations,
  plumes = SimulationsSmall@plumes[1:6,],
  values = raster(x = detectionMatrix2,
                  xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                  crs = "+init=epsg:4326"))

# test case: 2 locations are sufficient, but it probabily returns "4 are sufficient" and then "3 are sufficient"
detectionMatrix3 = matrix(c(1,1,1,1,0,0,0,0,
                            1,1,0,0,1,1,0,0,
                            1,0,1,0,1,0,1,0,
                            0,1,0,1,0,1,0,1),
                          nrow = 4, ncol = 8, byrow = TRUE)

completeExample3 = Simulations( 
  locations = SimulationsSmall@locations[1:4,],
  plumes = SimulationsSmall@plumes[c(1:5, 1:3),],
  values = raster(x = detectionMatrix3,
                  xmn = -90, xmx = 90, ymn = -90, ymx = 90,
                  crs = "+init=epsg:4326"))

# ------------ tests ---------------------------------

# test case from book: SD, cost correct
# no "nameSave" -> no files generated
filesBefore = list.files()
optSD_global1 = optimiseSD_global(
  simulations = completeExample, 
  aimNumber = 3,
  verbatim = FALSE,
  detectable = 1,          
  findAllOptima = TRUE     
)

expect_equal(# expected SD
  optSD_global1$SD,
  matrix(c(1, 2, 5,
           5, 6, 7),
         ncol = 3, byrow = TRUE)
)
expect_equal(# expected cost
  optSD_global1$cost,
  0
)
expect_equal(# no "nameSave" -> no files generated
  filesBefore, 
  list.files()
)


# non-complient data cause errors
# 'detect' selects layers
expect_error(
  optimiseSD_global(
  simulations = completeExample, 
  aimNumber = 3,
  verbatim = FALSE,
  detectable = 3,  # NA        
  findAllOptima = TRUE)
)  
expect_error(
  optimiseSD_global(
    simulations = completeExample, 
    aimNumber = 3,
    verbatim = FALSE,
    detectable = 4,  # != 0 1        
    findAllOptima = TRUE)
) 
expect_error(
  optimiseSD_global(
    simulations = completeExample, 
    aimNumber = 3,
    verbatim = FALSE,
    detectable = 5,  # NaN        
    findAllOptima = TRUE)
)

# logical or 0/1: no difference
# nameSave given -> files created (all there)
#filesToDelete =  list.files()[sapply(X = strsplit(list.files(), "optSD_global2"), FUN = "[[", 1) == ""]
#file.remove(filesToDelete)
filesBefore = list.files()
optSD_global2 = optimiseSD_global(
  simulations = completeExample, 
  aimNumber = 3,
  nameSave = "optSD_global2", 
  verbatim = FALSE,
  detectable = 2,          
  findAllOptima = TRUE     
)
newFiles = setdiff(list.files(), filesBefore)
expect_false(# there are new files
  identical(character(0), newFiles)  
)
expect_equal(# all new files are saved at 'saveName'
  sapply(X = strsplit(newFiles, "optSD_global2"), FUN = "[[", 1),
  rep("", 14)  
)
# saved data
load("optSD_global2_all_finalSD.Rdata")
finalSD
load("optSD_global2_first_1_SD_nRnC_1.Rdata")


expect_equal(# logical or 0/1: no difference 
  optSD_global1[1:2], # exclude report with Time (that differs)
  optSD_global2[1:2]
)
expect_equal(# logical or 0/1: no difference 
  optSD_global1[[3]][[2]][1:3], # exclude report with Time (that differs)
  optSD_global2[[3]][[2]][1:3]
)
expect_equal(# logical or 0/1: no difference 
  optSD_global1[[3]][[3]][[1]][1:3], # exclude report with Time (that differs)
  optSD_global2[[3]][[3]][[1]][1:3]
)
expect_equal(# logical or 0/1: no difference 
  optSD_global1[[3]][[3]][[2]][1:3], # exclude report with Time (that differs)
  optSD_global2[[3]][[3]][[2]][1:3]
)
expect_equal(# logical or 0/1: no difference 
  optSD_global1[[3]][[4]][1:3], # exclude report with Time (that differs)
  optSD_global2[[3]][[4]][1:3]
)

# reconstruction and back-transformation as expected
optSD_global3 = optimiseSD_global(
  simulations = completeExample2, 
  aimNumber = 3,
  nameSave = "optSD_global1", 
  verbatim = TRUE,
  detectable = "layer.2",          
  findAllOptima = TRUE    
)
expect_equivalent(
  optSD_global3$report$first$finalSD,
  matrix(c("4", "1", "3",
           "1", "3", "5"),
         ncol = 3, byrow = TRUE)
)
expect_equivalent(
  optSD_global3$report$all$finalSD,
  matrix(c("1", "3", "5",
           "5", "8", "9"),
         ncol = 3, byrow = TRUE)
)
expect_equal(
  optSD_global3$SD,
  rbind(
    matrix(c(rep(c(1,2), each = 2),
             rep(3,4),
             rep(c(5,7), times = 2)),
             ncol = 3),
    matrix(c(5,7,
             rep(8, 2),
             rep(9, 2)),
           ncol = 3)
  )
)
expect_equal(
  optSD_global3$cost,
  0
)


# also if some locations are deleted beforehand (locationsAll)
locationsAll2 = c(1,3,4,6,5,8,9) # turns completeExample2 into completeExample
optSD_global4 = optimiseSD_global(
  simulations = completeExample2, 
  aimNumber = 3,
  nameSave = "optSD_global4",
  verbatim = TRUE,
  detectable = "layer.2",          
  findAllOptima = TRUE,
  locationsAll = locationsAll2
)
expect_equal(
  optSD_global4$report$first$finalSD,
  optSD_global2$report$first$finalSD
)
expect_equal(
  optSD_global4$report$all$finalSD,
  optSD_global2$report$all$finalSD
)
expect_equal(
  optSD_global4$SD,
  matrix(c(1,3,5,
           5,8,9),
         nrow = 2, byrow = TRUE)
)
expect_equal(
  optSD_global4$cost,
  0
)


# locFix, locAll: deleting of rows and columns correct
expect_warning(
  optSD_global5 = optimiseSD_global(
    simulations = SimulationsSmall, 
    locationsAll = 1:8,
    locationsFix = 1,
    aimNumber = 1,
    verbatim = FALSE,
    detectable = "detect",          
    findAllOptima = FALSE     
  )  
)

expect_equivalent(
  optSD_global5$report$detectable,
  matrix(getValues(SimulationsSmall@values)[, "detect"], byrow = TRUE, nrow = 9)[2:8,c(2,5)]
  # plume 1 detected by fix sensor at location 1
  # plume 3 not detected at any location except 9 which is not in locationsAll
  # plume 4 not detected at all
)
expect_equal(
  optSD_global5$SD,
  matrix(4)        # reconstructed from locationsAll
)

#- - - - - - - locAll / locFix -> detectable has only 1 column or row - - - - - - -
# 1 column, aimNumber correct, all optima
optSD_global6a = optimiseSD_global(# 1 column (1st)
  simulations = SimulationsSmall, 
  locationsAll = 1:8,
  locationsFix = 4,
  aimNumber = 1,
  verbatim = FALSE,
  detectable = "detect",          
  findAllOptima = TRUE     
)
expect_equivalent(
  optSD_global6a$report$detectable,
  matrix(getValues(SimulationsSmall@values)[cellFromCol(SimulationsSmall@values, 1), "detect"], ncol = 1)[c(1:3, 5:8),,drop = FALSE]
)
expect_equal(# all optima
  optSD_global6a$SD,
  matrix(c(1,7), ncol = 1)
)
expect_equal(
  optSD_global6a$cost,
  0.4 # of 5 plumes 3 and 4 not detected (4 not detectable, 3 only at loc 9 which is not in locationsAll)  
)

# 1 column, aimNumber correct, one true optimum
optSD_global6b = optimiseSD_global(# 1 column (1st)
  simulations = SimulationsSmall, 
  locationsAll = 1:8,
  locationsFix = 4,
  aimNumber = 1,
  verbatim = FALSE,
  detectable = "detect",
  findSensorNumber = TRUE
)
expect_equal(# 1 optimum
  optSD_global6b$SD,
  matrix(1, ncol = 1) 
)
expect_equal(
  optSD_global6b$cost,
  0.4
)


# 1 column, aimNumber too high, all optima
optSD_global6c = optimiseSD_global(# 1 column
  simulations = SimulationsSmall, 
  locationsAll = 1:8,
  locationsFix = 4,
  aimNumber = 2, # more than dim of detectable
  verbatim = FALSE,
  detectable = "detect",          
  findAllOptima = TRUE     
)
expect_equal(# all optima
  optSD_global6c$SD,
  optSD_global6a$SD
)
expect_equal(
  optSD_global6c$cost,
  0.4
)

# 1 column, aimNumber too high, one true optimum
optSD_global6d = optimiseSD_global(# 1 column
  simulations = SimulationsSmall, 
  locationsAll = 1:8,
  locationsFix = 4,
  aimNumber = 2, # more than dim of detectable
  verbatim = FALSE,
  detectable = "detect",
  findSensorNumber = TRUE
)
expect_equal(# one optimum
  optSD_global6d$SD,
  optSD_global6b$SD
)
expect_equal(
  optSD_global6d$cost,
  0.4
)


# 1 column, aimNumber too high, one optimum
optSD_global6e = optimiseSD_global(# 1 column
  simulations = SimulationsSmall, 
  locationsAll = 1:8,
  locationsFix = 4,
  aimNumber = 2, # more than dim of detectable
  verbatim = FALSE,
  detectable = "detect"
)
expect_equal(# one optimum
  optSD_global6e$SD,
  optSD_global6b$SD
)
expect_equal(
  optSD_global6e$cost,
  0.4
)


# 1 location
# 1 optimum
optSD_global6f = optimiseSD_global(
  simulations = SimulationsSmall, 
  locationsAll = 4, # 1 location (with 2/5 detectable plumes)
  aimNumber = 1,
  verbatim = FALSE,
  detectable = "detect"    
)
expect_equal(
  optSD_global6f$SD,
  matrix(4)
)
expect_equal(
  optSD_global6f$cost,
  0.6
)  

# 1 optimum
optSD_global6g = optimiseSD_global(
  simulations = SimulationsSmall, 
  locationsAll = 3:6, # only 1 location with 2/5 detectable plumes
  aimNumber = 2,# more than needed
  verbatim = FALSE,
  detectable = "detect"    
)
expect_equal(
  optSD_global6g$SD,
  matrix(4) 
)
expect_equal(
  optSD_global6g$cost,
  0.6
)  


optSD_global6h = optimiseSD_global(
  simulations = SimulationsSmall, 
  locationsAll = 3:6, # only 1 location with detectable plumes
  aimNumber = 2,# more than needed
  verbatim = FALSE,
  detectable = "detect",
  findSensorNumber = TRUE
)
expect_equal(
  optSD_global6h$SD,
  matrix(4) 
)
expect_equal(
  optSD_global6h$cost,
  0.6
)  


expect_warning(
  expect_error(# no locations
    optimiseSD_global(
      simulations = SimulationsSmall,  
      locationsFix = 1:9, # no other locations left
      aimNumber = 1,
      verbatim = FALSE,
      detectable = "detect",          
      findAllOptima = TRUE     
    )
  )  
)

expect_warning(
  expect_error(# no detectable plume
    optimiseSD_global(
      simulations = SimulationsSmall, 
      locationsAll = c(3, 5, 6), # locations without detectable plumes
      aimNumber = 1,
      verbatim = FALSE,
      detectable = "detect"    
    )
  )
)
expect_warning(
  expect_error(# more sensors wanted than locations
    optimiseSD_global(
      simulations = SimulationsSmall, 
      locationsAll = 4, 
      aimNumber = 2, # more than locations
      verbatim = FALSE,
      detectable = "detect"    
    )
  )  
)
expect_warning(
  expect_error(# no locations
    optimiseSD_global(
      simulations = SimulationsSmall, 
      locationsAll = 10, # invalid -> no locations
      aimNumber = 2, # more than locations
      verbatim = FALSE,
      detectable = "detect"    
    )
  )
)

# --------------- less sensors are sufficient (but are not found in first search) ------------------
# all optima
optSD_global7a = optimiseSD_global(
  simulations = completeExample3, 
  aimNumber = 4,
  verbatim = TRUE,   
  findAllOptima = TRUE
)
expect_equivalent(
  optSD_global7a$SD,
  matrix(c(3,4), nrow = 1)
)  
expect_equal(
  optSD_global7a$cost,
  0
)
# one true optimum
optSD_global7b = optimiseSD_global(
  simulations = completeExample3, 
  aimNumber = 4,
  verbatim = TRUE,   
  findSensorNumber = TRUE
)
expect_equivalent(
  optSD_global7b$SD,
  matrix(c(3,4), nrow = 1)
)  
expect_equal(
  optSD_global7b$cost,
  0
)
# one "optimum" / one sensor set of desired length that detects all plumes
optSD_global7c = optimiseSD_global(
  simulations = completeExample3, 
  aimNumber = 4
)
expect_equivalent(
  optSD_global7c$SD,
  matrix(c(1,2,3,4), nrow = 1) # not a true optimum (as expected)!
)  
expect_equal(
  optSD_global7c$cost,
  0
)
# one sensor set of desired length that detects all plumes
optSD_global7d = optimiseSD_global(
  simulations = completeExample3, 
  aimNumber = 2
)
expect_equivalent(
  optSD_global7d$SD,
  matrix(c(3,4), nrow = 1) # true optimum (desired number of sensors fits)
)  
expect_equal(
  optSD_global7d$cost,
  0
)

# stop early: maxIterations
optSD_global8a = optimiseSD_global(
  simulations = completeExample, 
  aimNumber = 3,
  detectable = 1,          
  findAllOptima = TRUE,
  maxIterations = choose(7,3) # all possibilities
)
expect_equivalent(
  optSD_global8a$SD,
  optSD_global1$SD
)  

optSD_global8b = optimiseSD_global(
  simulations = completeExample, 
  aimNumber = 3,
  detectable = 1,          
  findAllOptima = TRUE,
  maxIterations = 5 # less than all possibilities
)
expect_equivalent(
  optSD_global8b$SD, # solutions that find 6 plumes (limit derived from first run)
  matrix(c(rep(3, 2 * 3 +          1 * 2),
           rep(c(1,7), each = 3),  c(1,7),
         rep(c(2,4,6), times = 2), rep(5, 2)),
         ncol = 3)
) 

optSD_global8c = optimiseSD_global(
  simulations = completeExample, 
  aimNumber = 3,
  detectable = 1,          
  findSensorNumber = TRUE,
  maxIterations = 5 # less than all possibilities
)
expect_equivalent(
  optSD_global8c$SD,
  optSD_global8b$SD[1,,drop = FALSE]
)  

# get aimNumber from locInit (after correction)
optSD_global9a = optimiseSD_global(
  simulations = SimulationsSmall, 
  locationsAll = 1:8,
  locationsFix = 1,
  locationsInitial = 1:2,
  detectable = "detect",          
  findAllOptima = FALSE     
)  
expect_equal(
  optSD_global9a[1:2],
  optSD_global5[1:2]
)
expect_equal(
  optSD_global9a[[3]][[2]][[2]],
  optSD_global5[[3]][[2]][[2]]
)

expect_error(
  optimiseSD_global(# no aimNumber, no valid locationsInitial -> aimNumber = 0 -> error
    simulations = SimulationsSmall, 
    locationsAll = 1:8,
    locationsFix = 1,
    locationsInitial = 1,
    detectable = "detect",          
    findAllOptima = FALSE     
  )    
)

# example of realistic size
optSD_global10 = optimiseSD_global(
  simulations = radioactivePlumes_area, 
  locationsAll = locAll_l,
  aimNumber = 3,
  detectable = "detectable",          
  findAllOptima = TRUE,
  nameSave = "optSD_global9" 
)
SD10 = as(subsetSDF(radioactivePlumes_area@locations,
                   locations =  optSD_global10$SD[1,]),
         "SpatialPointsDataFrame")
spplot(subsetSDF(radioactivePlumes_area@locations, locations = locAll_l),
       sp.layout = list("sp.points", 
                        SpatialPoints(coords = coordinates(SD10), proj4string = CRS(proj4string(SD10))),
                        col = 3, cex = 2, lwd = 1.5)
)

# ---------- call via optimiseSD ------------------
optSDglobal = replaceDefault(
  fun = optimiseSD_global,
  newDefaults = list(
    detectable = 1,          
    findAllOptima = TRUE 
    ),
  type = "optimisationFun.optimiseSD"
  )

optSD_global1 = optimiseSD_global(
  simulations = completeExample, 
  aimNumber = 3,
  verbatim = FALSE,
  detectable = 1,          
  findAllOptima = TRUE     
)

optSD_global11 = optimiseSD( 
  simulations = completeExample,
  aimNumber = 3,
  costFun = NA,
  optimisationFun = optSDglobal[[1]],
  nameSave = NA
)

expect_equal(
  optSD_global1[1:2],
  optSD_global11[1:2]
)


# completeExample4
#  matrix(c(rep(c(1,2,3,9), each = 3 * 2), # detect 1,2,5,6
#           rep(rep(c(4,6,8), each = 2), times = 4), # detect 1,2,3,4
#           rep(c(5,7), times = 4 * 3)), # detect 1,3,5,7
#  ncol = 3)
