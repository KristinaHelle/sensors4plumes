##################################################
# test subset.Simulations                        #
##################################################

test_that("subset.Simulations location order", {
  data(radioactivePlumes_area)
  data(radioactivePlumes_local)
  # case: input and output is in memory
  ## no subsetting
  simulations_memory_all = subset(radioactivePlumes_local)
  expect_equal(
    simulations_memory_all,
    radioactivePlumes_local
  )
  ### valuesOnly
  simulations_memory_all_vO = subset(radioactivePlumes_local, valuesOnly = TRUE)
  expect_equal(
    simulations_memory_all_vO,
    radioactivePlumes_local@values
  )  
  ## sorted, small indices
  locations = sort(sample.int(250, 5))
  plumes = sort(sample.int(200, 5))
  simulations_memory = subset(radioactivePlumes_local, locations = 1:250, plumes = 1:200)

  subset_memory_l = subset(simulations_memory, locations = locations)
  subset_memory_p = subset(simulations_memory, plumes = plumes)
  subset_memory_pl = subset(simulations_memory, plumes = plumes, locations = locations) 
  
  expect_equal(
    getValues(subset_memory_l@values, 3, 1),
    simulations_memory@values[locations[3],]
  )     
  expect_equal(
    getValues(subset_memory_p@values, locations[1], 1),
    simulations_memory@values[locations[1], plumes]
  ) 
  expect_equal(
    getValues(subset_memory_pl@values, 1, 1),
    simulations_memory@values[locations[1],plumes]
  ) 
  
  ## indices with switched order
  locationsR = sample.int(nLocations(radioactivePlumes_local), 500)
  plumesR = sample.int(nPlumes(radioactivePlumes_local), 400)
  subset_memory_R = subset(radioactivePlumes_local, locations = locationsR, plumes = plumesR)
  expect_equal(
    getValues(subset_memory_R@values, locations[4], 1),
    getValues(radioactivePlumes_local@values, locationsR[locations[4]], 1)[plumesR,]
  )
  ## indices with invalid values
  ### normal mode: invalid indices ignored
  subset_memory_RNA = subset(radioactivePlumes_local, 
                               locations = c(0, NA, 10000, 0, locationsR), 
                                             plumes = c(plumesR, 0, NA, 10000, 0))
  expect_equal(
    subset_memory_R,
    subset_memory_RNA
  )  
  
  locationsR2 = sample.int(nLocations(simulations_memory), 250)
  plumesR2 = sample.int(nPlumes(simulations_memory), 200)
  subset_memory_RNA2 = subset(radioactivePlumes_local, 
                                 locations = c(0, NA, 10000, 0, locationsR2), 
                                 plumes = c(plumesR2, 0, NA, 10000, 0))
  
  ### valuesOnly mode: invalid indices replaced by NA values
  subset_memory_RNA_vO = subset(simulations_memory, 
                                 locations = c(0, NA, 10000, 0, locationsR2), 
                                 plumes = c(plumesR2, 0, NA, 10000, 0),
                                 valuesOnly = TRUE) 
  expect_true(
    is(subset_memory_RNA_vO, "RasterStack"))
  expect_equivalent(
    subset_memory_RNA_vO[1:4,],
    matrix(NA, nrow = 4 * 204, ncol = 3)
  )
  expect_equivalent(
    subset_memory_RNA_vO[,201:204],
    matrix(NA, nrow = 4 * 254, ncol = 3)
  )
  expect_equivalent(
    subset_memory_RNA_vO[5:254, 1:200],
    getValues(subset_memory_RNA2@values)
  )
  
  ## indices with repetitions
  ### normal mode: repetitions ignored
  subset_memory_R2 = subset(simulations_memory, 
                              plumes = c(plumes, plumes), 
                              locations = c(locations, locations))
  expect_equal(
    subset_memory_R2,
    subset_memory_pl
  )   
  ### valuesOnly mode: repetitions taken into account
  subset_R_vO = subset(simulations_memory, 
                       locations = c(locations, locations), 
                       plumes = c(plumes, plumes), 
                       valuesOnly = TRUE) 
  expect_true(
    is(subset_R_vO, "RasterStack"))
  expect_equal(
    subset_R_vO[1:5,1:5],
    getValues(subset_memory_pl@values)
  )
  expect_equal(
    subset_R_vO[6:10,1:5],
    getValues(subset_memory_pl@values)
  )
  expect_equal(
    subset_R_vO[1:5,6:10],
    getValues(subset_memory_pl@values)
  )
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - 
  # case: input is processed in chunks, but not output
  ## no subsetting
  simulations_memory_all2 = subset(radioactivePlumes_area)
  expect_equal(
    simulations_memory_all2,
    radioactivePlumes_area
  )
  ### valuesOnly
  simulations_memory_all2_vO = subset(radioactivePlumes_area, valuesOnly = TRUE)
  expect_equal(
    simulations_memory_all2_vO,
    radioactivePlumes_area@values
  ) 
  ### normal case, multiple indices ignored
  locations2 = sample(nLocations(radioactivePlumes_area), 15)
  plumes2 = sample(nPlumes(radioactivePlumes_area), 15)
  i = sample.int(15,1)
  subset_fm_RNA = subset(radioactivePlumes_area, 
                        plumes = c(plumes2, NA, 10000, NA, 0, plumes2), 
                        locations = c(locations2, NA, 10000, NA, 0,locations2))
  expect_true(
    inMemory(subset_fm_RNA@values)
  )
  expect_equal(
    getValues(subset_fm_RNA@values, i, 1),
    getValues(radioactivePlumes_area@values, locations2[i], 1)[plumes2, ]
  )
  expect_equal(
    dim(subset_fm_RNA@values),
    c(length(locations2), length(plumes2), 3)
  )
  
  ### valuesOnly: invalid indices -> repetitions taken into account
  result_mf_R_vO = subset(radioactivePlumes_area, 
                         plumes = c(plumes2, plumes2), 
                         locations = c(locations2, locations2),
                         valuesOnly = TRUE)
  expect_equal(
    getValues(result_mf_R_vO, i + 15, 1),
    rbind(getValues(subset_fm_RNA@values, i, 1), 
          getValues(subset_fm_RNA@values, i, 1))
  )
  expect_equal(
    dim(result_mf_R_vO),
    c(2 * length(locations2), 2 * length(plumes2), 3)
  )
  ### valuesOnly: invalid indices -> NA taken into account
  subset_mf_RNA_vO = subset(radioactivePlumes_area, 
                             plumes = c(plumes2, NA, 0, 10000, NA, plumes2), 
                             locations = c(locations2, NA, NA, 10000, 0, 0, locations2),
                             valuesOnly = TRUE)
  expect_equal(
    getValues(subset_mf_RNA_vO, i + 15 + 5, 1),
    rbind(getValues(subset_fm_RNA@values, i, 1), NA, NA, NA, NA, 
          getValues(subset_fm_RNA@values, i, 1))
  )
  expect_equivalent(
    getValues(subset_mf_RNA_vO, 16, 4),
    matrix(as.numeric(NA), nrow = 4 * 34, ncol = 3)
  )
  
  # case: input and output is processed in chunks, output is saved to disk
  ## all data, but changed order
  ### normal case: order taken into account
  random_locationsAll = sample(nLocations(radioactivePlumes_area))
  random_plumesAll = sample(nPlumes(radioactivePlumes_area))
  expect_warning(
    subset_file <- subset(radioactivePlumes_area, 
                         plumes = random_plumesAll, 
                         locations = random_locationsAll,
                         overwrite = TRUE))
  
  expect_false(inMemory(subset_file@values))
  
  i = sample.int(2500, 1)
  expect_equal(
    getValues(subset_file@values, i, 1),
    getValues(radioactivePlumes_area@values, random_locationsAll[i], 1)[random_plumesAll, ]
  )
  ## with repetitions
  ### normal case: NA, repetitions ignored
  random_locationsMult = sample(nLocations(radioactivePlumes_area), 2000, replace = TRUE)
  random_plumesMult = sample(nPlumes(radioactivePlumes_area), 2000, replace = TRUE)
  expect_warning(
    subset_file_RNA <- subset(radioactivePlumes_area, 
                          plumes = c(NA, 0, 10000, random_plumesMult), 
                          locations = c(10000, 0, 0, random_locationsMult, NA),
                          overwrite = TRUE)
  )
  
  i = sample.int(length(unique(random_locationsMult)), 1)
  expect_equal(
    getValues(subset_file_RNA@values, i, 1),
    getValues(radioactivePlumes_area@values, 
              unique(random_locationsMult)[i], 1)[unique(random_plumesMult), ])  
  
  ### valuesOnly: repetition taken into account
  expect_warning(
    subset_file_R_vO <- subset(radioactivePlumes_area, 
                              plumes = random_plumesMult, 
                              locations = random_locationsMult,
                              valuesOnly = TRUE,
                              overwrite = TRUE)
    )
  expect_equal(
    getValues(subset_file_R_vO, i, 1),
    getValues(radioactivePlumes_area@values, 
              random_locationsMult[i], 1)[random_plumesMult, ])  
 
  ### valuesOnly: repetition taken into account
  expect_warning(
    result_file_RNA_vO <- subset(radioactivePlumes_area, 
                                 plumes = c(NA, 0, 10000, random_plumesMult), 
                                 locations = random_locationsMult,
                                 valuesOnly = TRUE,
                                 overwrite = TRUE)
  )
  
  expect_equal(
    getValues(result_file_RNA_vO, i, 1),
    getValues(radioactivePlumes_area@values, 
              random_locationsMult[i], 1)[c(NA, NA, NA, random_plumesMult), ])  
  expect_equal(
    dim(result_file_RNA_vO),
    c(2000, 2003, 3)
  )
  
  expect_warning(
    result_file_RNA_vO2 <- subset(radioactivePlumes_area, 
                                   plumes = c(NA, 0, 10000, random_plumesMult), 
                                   locations = c(10000, 0, 0, random_locationsMult, NA, NA),
                                   valuesOnly = TRUE,
                                   overwrite = TRUE)    
  )
  expect_equal(
    dim(result_file_RNA_vO2),
    c(2005, 2003, 3)
  )
  expect_equivalent(
    getValues(result_file_RNA_vO2, row = 1, nrows = 3),
    matrix(NA, nrow = 3 * 2003, ncol = 3)
  )  
  expect_equivalent(
    getValues(result_file_RNA_vO2, row = 2004, nrows = 2),
    matrix(NA, nrow = 2 * 2003, ncol = 3)
  )  
  
  expect_equal(
    getValues(result_file_RNA_vO2, i, 1),
    getValues(radioactivePlumes_area@values, 
              random_locationsMult[i-3], 1)[c(NA, NA, NA, random_plumesMult), ]) # i-3 to account for NA lines
  
  # case empty locations or plumes
  expect_error(
    subset(radioactivePlumes_area, plumes = c(5000, 6000))
  )
  expect_error(
    subset(radioactivePlumes_area, locations = c(7000, 8000))
  )
)

