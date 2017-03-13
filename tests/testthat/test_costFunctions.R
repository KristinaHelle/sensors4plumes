########################################################################
# test cost functions                                                  #
########################################################################
# turn all pre cost functions into real cost functions and test them
## interpolationError
# absError
# delineationError
# absErrorMap
# delineationErrorMap
## minimalDistance
# krigingVariance
# spatialSpread
## measurementsResult
# singleDetection
# multipleDetection
# earlyDetection


data(SimulationsSmall)
locDel1 = sample.int(nLocations(SimulationsSmall), 4)
locKeep1 = sample(setdiff(1:nLocations(SimulationsSmall), locDel1), 2)
data(radioactivePlumes_local)
locDel2 = sample.int(nLocations(radioactivePlumes_local), 5)
locKeep2 = sample(setdiff(1:nLocations(radioactivePlumes_local), locDel2), 100)
data(radioactivePlumes_area)
locDel3 = sample.int(nLocations(radioactivePlumes_area), 5)
locKeep3 = sample(setdiff(1:nLocations(radioactivePlumes_area), locDel3), 100)

#---------------------------------------------------------------------------#
                        # interpolationError
#---------------------------------------------------------------------------#
### interpolation functions
#### median variogram
##### SimulationsSmall: error (autofitVariogram): needs square grid 
# medianVariogram_s1 = fitMedianVariogram(simulations = SimulationsSmall, 
#                                         plumes = 1:nPlumes(SimulationsSmall),
#                                         values = 1) 
##### radioactivePlumes_local
medianVariogram_l1 = fitMedianVariogram(simulations = radioactivePlumes_local, 
                                         plumes = sample.int(nPlumes(radioactivePlumes_local), 100),
                                         values = 1) 
##### radioactivePlumes_area
medianVariogram_a1 = fitMedianVariogram(simulations = radioactivePlumes_area, 
                                        plumes = sample.int(nPlumes(radioactivePlumes_area), 100),
                                        values = 1)
save(medianVariogram_l1, medianVariogram_a1, file = "/home/kristina/Desktop/s4p_extra/data/medianVariograms.Rdata")

#### interpolation functions
##### kriging radioactivePlumes_local
krige0var_l1 = replaceDefault(krige0, newDefaults = list(
  formula = z ~ 1, model = medianVariogram_l1, beta = NA, ... = NA))[[1]]
##### kriging radioactivePlumes_area
krige0var_a1 = replaceDefault(krige0, newDefaults = list(
  formula = z ~ 1, model = medianVariogram_a1, beta = NA, ... = NA))[[1]]
##### inverse distance weighting
idw0z = replaceDefault(idw0, newDefaults = list(
  formula = z ~ 1))[[1]]

# - - -- - - - -- - - delination error -- - - - - - -- - - #
# define cost function
delineationError_kl1 = replaceDefault(
  interpolationError, newDefaults = list(
    values = 1,
    fun_interpolation = krige0var_l1,
    fun_error = delineationError, # must return only 1 value
    fun_Rpl_cellStats = "mean",
    tmpfile = "interpolationError_del_kl1"), 
  type = "costFun.optimiseSD")[["fun"]]

# compute cost
delineationError_kl1_1 = 
  delineationError_kl1(simulations = radioactivePlumes_local,
                       locations = c(locDel2[-1], locKeep2))
# apply cost function in 'deleteSensor'
del_int_kl1 = deleteSensor(simulations = radioactivePlumes_local, 
                    costFun = delineationError_kl1, 
                    locationsDeletable = locDel2, 
                    locationsKeep = locKeep2)
# check result of 'deleteSensor'
costWithout_del_int_kl1 = numeric(length(locDel2))
for (i in seq(along = locDel2)){
  costWithout_del_int_kl1[i] = delineationError_kl1(radioactivePlumes_local, c(locDel2[-i], locKeep2))[[1]]
}
expect_equal(
  del_int_kl1[[3]],
  costWithout_del_int_kl1)
expect_equal(
  del_int_kl1[[2]],
  which(costWithout_del_int_kl1 == min(costWithout_del_int_kl1))
)
expect_true(
  is.element(setdiff(locDel2, del_int_kl1[[1]]), 
             locDel2[which(costWithout_del_int_kl1 == min(costWithout_del_int_kl1))])
)

# delineationErrorMap
# -  -- - - - - - - - - -- absolute error  -  - --  -- -- - - - - -- #
absError_ka1 = replaceDefault(
  interpolationError, newDefaults = list(
    values = 1,
    fun_interpolation = krige0var_a1,
    fun_error = absError, 
    fun_Rpl_cellStats = "mean",
    tmpfile = "interpolationError_abs_ka1"
  ), 
  type = "costFun.optimiseSD")[["fun"]]

absError_i1 = replaceDefault(
  interpolationError, newDefaults = list(
    values = 1,
    fun_interpolation = idw0z,
    fun_error = absError, # must return only 1 value
    fun_Rpl_cellStats = "mean",
    tmpfile = "interpolationError_del_i1"), 
  type = "costFun.optimiseSD")[["fun"]]


# compute cost
absError_ka1_1 = 
  absError_ka1(simulations = radioactivePlumes_area,
                       locations = c(locDel3[-1], locKeep3))
# apply cost function in 'deleteSensor'
del_int_kl1 = deleteSensor(simulations = radioactivePlumes_area, 
                           costFun = absError_ka1, 
                           locationsDeletable = locDel3, 
                           locationsKeep = locKeep3)


# absErrorMap

# --------------------------------------------------------------- #
                    # spatial spread 
#-----------------------------------------------------------------#
# - - -  - - -- - - - minimal distance - - - - - - ---- - - - - -- #
# define cost function
meanFun = function(x){mean(x, na.rm = TRUE)}
minDist = replaceDefault(
  spatialSpread, newDefaults = list(
    fun = minimalDistance,
    fun_R = meanFun
  ), type = "costFun.optimiseSD"
)[["fun"]] 

# compute cost
minDist1 = minDist(simulations = SimulationsSmall,
                   locations = c(locDel1[-1], locKeep1))
# apply cost function in 'deleteSensor'
del_minDist_1 = deleteSensor(simulations = SimulationsSmall, 
                    costFun = minDist, 
                    locationsDeletable = locDel1, 
                    locationsKeep = locKeep1)


# - - --  -- - -- - - - -krigingVariance - - - - - --- - - - - -#
# define cost function
krigVar_a1 = replaceDefault(
  spatialSpread, newDefaults = list(
    fun = replaceDefault(krigingVariance, newDefaults = list(model = medianVariogram_a1))[["fun"]],
    fun_R = meanFun
  ), type = "costFun.optimiseSD"
)[["fun"]] 

# compute cost
krigVar_a1_1 = krigVar_a1(simulations = radioactivePlumes_area,
                   locations = c(locDel3[-1], locKeep3))

# apply cost function in 'deleteSensor'
del_krigVar_a1_1 = deleteSensor(simulations = radioactivePlumes_area, 
                             costFun = krigVar_a1, 
                             locationsDeletable = locDel3, 
                             locationsKeep = locKeep3)

# ----------------------------------------------------------------------- #
                        # measurementsResult
# ----------------------------------------------------------------------- #
# generate intermediate results
## layer: detection 
threshold = 1e-7
radioactivePlumes_local@values$detectable = calc(
  radioactivePlumes_local@values$maxdose,
  fun = function(x){x >= threshold})
## plumes: 
### total dose
radioactivePlumes_local@plumes$totalDose = 
  summaryPlumes(radioactivePlumes_local, fun = sum, values = "finaldose")[[2]]
### number of locations where it can be detected
radioactivePlumes_local@plumes$nDetectable = 
  summaryPlumes(radioactivePlumes_local, fun = sum, values = "detectable")[[2]]  
### earliest possible detection of plume inside the area
radioactivePlumes_local@plumes$earliestDetection = 
  

#  - - - - -- - - - - - - singleDetection - - - - - - -- - - - - -- - --  #
# define cost function: see measurementsResult.R
# compute cost
sDet_1 = singleDetection(simulations = radioactivePlumes_local,
                         locations = c(locDel2[-1], locKeep2))
# apply cost function in 'deleteSensor'
del_sDet_1 = deleteSensor(simulations = radioactivePlumes_local, 
                                costFun = singleDetection, 
                                locationsDeletable = locDel2, 
                                locationsKeep = locKeep2)
# check result of 'deleteSensor'
costWithout_sDet_1 = numeric(length(locDel2))
for (i in seq(along = locDel2)){
  costWithout_sDet_1[i] = singleDetection(radioactivePlumes_local, c(locDel2[-i], locKeep2))[[1]]
}
expect_equal(
  del_sDet_1[[3]],
  costWithout_sDet_1)
expect_equal(
  del_sDet_1[[2]],
  which(costWithout_sDet_1 == min(costWithout_sDet_1))
)
expect_true(
  is.element(setdiff(locDel2, del_sDet_1[[1]]), 
             locDel2[which(costWithout_sDet_1 == min(costWithout_sDet_1))])
)


#  - - - -- - - -- - multipleDetection - - - -- - -- ---- - - - - -- - - - #
# define cost function: see measurementsResult.R
# compute cost
mDet_1 = multipleDetection(simulations = radioactivePlumes_local,
                         locations = c(locDel2[-1], locKeep2))
# apply cost function in 'deleteSensor'
del_mDet_1 = deleteSensor(simulations = radioactivePlumes_local, 
                          costFun = multipleDetection, 
                          locationsDeletable = locDel2, 
                          locationsKeep = locKeep2)

# check result of 'deleteSensor'
costWithout_mDet_1 = numeric(length(locDel2))
for (i in seq(along = locDel2)){
  costWithout_mDet_1[i] = multipleDetection(radioactivePlumes_local, c(locDel2[-i], locKeep2))[[1]]
}
expect_equal(
  del_mDet_1[[3]],
  costWithout_mDet_1)
expect_equal(
  del_mDet_1[[2]],
  which(costWithout_mDet_1 == min(costWithout_mDet_1))
)
expect_true(
  is.element(setdiff(locDel2, del_mDet_1[[1]]), 
             locDel2[which(costWithout_mDet_1 == min(costWithout_mDet_1))])
)


# - - - - - - - - - --  earlyDetection - - - - -  - -- --- -- --  - #
# define cost function: see measurementsResult.R
# compute cost
eDet_1 = earlyDetection(simulations = radioactivePlumes_local,
                           locations = c(locDel2[-1], locKeep2))
# apply cost function in 'deleteSensor'
del_eDet_1 = deleteSensor(simulations = radioactivePlumes_local, 
                          costFun = earlyDetection, 
                          locationsDeletable = locDel2, 
                          locationsKeep = locKeep2)
