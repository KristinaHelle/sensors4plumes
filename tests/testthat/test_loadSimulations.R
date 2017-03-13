#########################################################################
# test                        loadSimulations                           #
#########################################################################

# A) raster files, a) per file: single plume, single kind
#basicPath1 = "/home/kristina/Desktop/s4p/data/fileFormats/test_raster"
basicPath1 = "/home/kristina/Desktop/s4p/data/fileFormats/raster"
files1 = list.files(basicPath1)

loadedSimulations1 = loadSimulations(
  basicPath = basicPath1,
  savePath = "/home/kristina/Desktop/s4p_article/test/t1_",
  overwrite = TRUE
)
test_that("loadSimulations raster kind in files",{
  # correct dimensions
  expect_equal(nLocations(loadedSimulations1), 100)
  expect_equal(nPlumes(loadedSimulations1), 5)
  expect_equal(nKinds(loadedSimulations1), 3)
  # correct values
  expect_equal(
    getValues(raster(paste0(basicPath1, "/", files1[1]))),
    loadedSimulations1@values[[1]][,1]
  )  
  expect_equal(
    getValues(raster(paste0(basicPath1, "/", files1[2]))),
    loadedSimulations1@values[[2]][,1]
  )  
  expect_equal(
    getValues(raster(paste0(basicPath1, "/", files1[4]))),
    loadedSimulations1@values[[1]][,2]
  ) 
})

# A)               b) single plume, all kinds
#basicPath2 = "/home/kristina/Desktop/s4p/data/fileFormats/test_rasterMultilayer"
basicPath2 = "/home/kristina/Desktop/s4p/data/fileFormats/raster_multilayer"
loadedSimulations2 = loadSimulations(
  basicPath = basicPath2,
  savePath = "/home/kristina/Desktop/s4p_article/test/t2_",
  overwrite = TRUE
)

test_that("loadSimulations raster kind multilayer",{
  # correct values (including correctly used names)
  expect_equal(dename(getValues(loadedSimulations1@values)[,c("typea", "typeb","typec")], kind = "dimnames"), 
               dename(getValues(loadedSimulations2@values)[,c("sourceA.date1_typea", "sourceA.date1_typeb", "sourceA.date1_typec")], kind = "dimnames")
  )  
})

# files with different number of layers
# basicPath2a = "/home/kristina/Desktop/s4p/data/fileFormats/raster_multilayer_plumes"
basicPath2a = "/home/kristina/Desktop/s4p/data/fileFormats/raster_multilayer_plumes"
 expect_error(
   loadedSimulations2a <- loadSimulations(
   basicPath = basicPath2a,
   savePath = "/home/kristina/Desktop/s4p_article/test/t2a_",
   overwrite = TRUE)
 )  


# A)               c) several plumes               
loadedSimulations3 = loadSimulations(# interpret kinds as plumes
  basicPath = basicPath2,
  savePath = "/home/kristina/Desktop/s4p_article/test/t3_",
  multilayer = "plumes",
  overwrite = TRUE
)               
test_that("loadSimulations raster several plumes in files",{
  # correct dimensions
  expect_equal(nLocations(loadedSimulations3), 100)
  expect_equal(nPlumes(loadedSimulations3), 15)
  expect_equal(nKinds(loadedSimulations3), 1)     
  # correct values (layer order as in loadedSimulations2)
  i = sample.int(5,1); j = sample.int(3,1)
  expect_equal(loadedSimulations3@values[,(i-1) * 3 + j][,1],
               loadedSimulations2@values[,i][,j])
})

# files with different number of layers
 loadedSimulations3a <- loadSimulations(
     basicPath = basicPath2a,
     savePath = "/home/kristina/Desktop/s4p_article/test/t3a_",
     multilayer = "plumes",
     overwrite = TRUE)
 test_that("loadSimulations raster several plumes in files of differing layer number",{
   # correct dimension
   expect_equal(nPlumes(loadedSimulations3a), 5)
   # correct values (layers shifted)
   expect_equal(loadedSimulations3a@values[,1][,1],
                loadedSimulations3@values[,3][,1])      
 })
                         
# A)               a)                           i) region
region = as.logical(replicate(100, sample.int(2,1)) - 1)
loadedSimulations4 = loadSimulations(
  basicPath = basicPath1,
  savePath = "/home/kristina/Desktop/s4p_article/test/t4_",
  overwrite = TRUE,
  region = region
)
test_that("loadSimulations raster kind in files, region",{
  # correct dimensions
  expect_equal(nLocations(loadedSimulations4), sum(region))
  expect_equal(nPlumes(loadedSimulations4), 5)
  expect_equal(nKinds(loadedSimulations4), 3)     
  # correct values
  i = sample.int(5,1); j = sample.int(3,1)
  expect_equal(loadedSimulations4@values[,i][,j],
               loadedSimulations1@values[,i][,j][region])
})

# A)               b)                           i) 
loadedSimulations5 = loadSimulations(
  basicPath = basicPath2,
  savePath = "/home/kristina/Desktop/s4p_article/test/t5_",
  overwrite = TRUE,
  region = region
)
test_that("loadSimulations raster kind multilayer, region",{
  expect_equal(nLocations(loadedSimulations5), sum(region))
  expect_equal(dename(getValues(loadedSimulations4@values), kind = "dimnames"), 
               dename(getValues(loadedSimulations5@values)[,c(3,2,1)], kind = "dimnames")
  )
})  

# A)               c)                           i) 
loadedSimulations6 = loadSimulations(# interpret kinds as plumes
  basicPath = basicPath2,
  savePath = "/home/kristina/Desktop/s4p_article/test/t6_",
  multilayer = "plumes",
  overwrite = TRUE,
  region = region
)               
test_that("loadSimulations raster several plumes in files, region",{
  # correct dimensions
  expect_equal(nLocations(loadedSimulations6), sum(region))
  expect_equal(nPlumes(loadedSimulations6), 15)
  expect_equal(nKinds(loadedSimulations6), 1)     
  # correct values
  i = sample.int(5,1); j = sample.int(3,1)
  expect_equal(loadedSimulations6@values[,(i-1) * 3 + j][,1],
               loadedSimulations5@values[,i][,j])
})

# A)               a)                                                  I) bbox
bboxOrig = bbox(loadedSimulations1@locations)
l = sort(sample.int(min(loadedSimulations1@locations@grid@cells.dim) - 1,4))
bBox = c(bboxOrig[1,1] + loadedSimulations1@locations@grid@cellsize[1] * l[1],
         bboxOrig[1,1] + loadedSimulations1@locations@grid@cellsize[1] * l[3],
         bboxOrig[2,1] + loadedSimulations1@locations@grid@cellsize[2] * l[2],
         bboxOrig[2,1] + loadedSimulations1@locations@grid@cellsize[2] * l[4])

loadedSimulations7 = loadSimulations(
  basicPath = basicPath1,
  savePath = "/home/kristina/Desktop/s4p_article/test/t7_",
  overwrite = TRUE,
  bBox = bBox
)

test_that("loadSimulations raster kind in files, bBox",{
  # correct dimensions
  expect_equal(nLocations(loadedSimulations7), (l[3] - l[1]) * (l[4] - l[2]))
  expect_equal(nPlumes(loadedSimulations7), 5)
  expect_equal(nKinds(loadedSimulations7), 3)
  # correct values
  i = sample.int(5,1); j = sample.int(3,1)
  expect_equal(
    getValues(crop(raster(paste0(basicPath1, "/", files1[(i-1) * 3 + j])), bBox)),
    loadedSimulations7@values[,i][,j]
  )
})

# A)               a)                           i)                     I) 
regionBBox = as.logical(matrix(region, nrow = 10)[(l[2] + 1):l[4] , (l[1] + 1):l[3]])
loadedSimulations8 = loadSimulations(
  basicPath = basicPath1,
  savePath = "/home/kristina/Desktop/s4p_article/test/t8_",
  overwrite = TRUE,
  region = regionBBox,
  bBox = bBox
)
# wrong region -> region ignored
expect_warning(
  loadedSimulations9 <- loadSimulations(
  basicPath = basicPath1,
  savePath = "/home/kristina/Desktop/s4p_article/test/t9_",
  overwrite = TRUE,
  region = c(regionBBox, TRUE),
  bBox = bBox)
)

test_that("loadSimulations raster kind in files, region, bBox",{
  # correct dimensions
  expect_equal(nLocations(loadedSimulations8), sum(regionBBox))
  expect_equal(nPlumes(loadedSimulations8), 5)
  expect_equal(nKinds(loadedSimulations8), 3)     
  # correct values
  i = sample.int(5,1); j = sample.int(3,1)
  
  expect_equal(
    getValues(crop(raster(paste0(basicPath1, "/", files1[(i-1) * 3 + j])), bBox))[regionBBox],
        loadedSimulations8@values[,i][,j]
    )
  expect_equal(getValues(loadedSimulations9@values),
               getValues(loadedSimulations7@values)
    )
})
                         
                         
# B) text files    a) 
basicPath3 = "/home/kristina/Desktop/s4p/data/fileFormats/text"
files3 = list.files(basicPath3)

LoadedSimulations10 = loadSimulations(
  basicPath = basicPath3,
  savePath = "/home/kristina/Desktop/s4p_article/test/t10_",
  overwrite = TRUE
)
loadedSimulations10 = Simulations(
  values = LoadedSimulations10[["values"]],
  plumes = LoadedSimulations10[["plumes"]],
  locations = loadedSimulations1@locations)

test_that("loadSimulations txt kind in files",{
  # correct dimensions
  expect_equal(nLocations(loadedSimulations10), 100)
  expect_equal(nPlumes(loadedSimulations1), 5)
  expect_equal(nKinds(loadedSimulations1), 3)
  # correct values
  expect_equal(
    loadedSimulations1@values[,1][,1],
    loadedSimulations10@values[,1][,1]
  )  
})

# B)               b)
basicPath4 = "/home/kristina/Desktop/s4p/data/fileFormats/text_multicolumn"
files4 = list.files(basicPath4)

LoadedSimulations11 = loadSimulations(
  basicPath = basicPath4,
  savePath = "/home/kristina/Desktop/s4p_article/test/t11_",
  overwrite = TRUE,
  readBy = "scan",
  skip = 1
) 
loadedSimulations11 = Simulations(
  values = LoadedSimulations11[["values"]],
  plumes = LoadedSimulations11[["plumes"]],
  locations = loadedSimulations1@locations)

test_that("loadSimulations txt kind multilayer",{
  expect_equal(dename(getValues(loadedSimulations11@values), kind = "dimnames"),
               dename(getValues(loadedSimulations10@values), kind = "dimnames")
  )
})  
# differnt header form
# basicPath5 = "/home/kristina/Desktop/s4p/data/fileFormats/test_txtMulticolumn_2lineHeader"
# files5 = list.files(basicPath5)
# LoadedSimulations12 = loadSimulations(
#   basicPath = basicPath5,
#   savePath = "/home/kristina/Desktop/s4p_article/test/t12_",
#   overwrite = TRUE,
#   readBy = "scan",
#   skip = 2
# ) 
# 
# test_that("loadSimulations txt kind multilayer 2line header",{
#   expect_equal(names(LoadedSimulations12[["values"]]),
#                names(read.table(paste0(basicPath5, "/", files5[1]), header = TRUE))
#   )        
#   expect_equal(nrow(LoadedSimulations12[["values"]]), 7629)
#   expect_equal(matrix(scan(paste0(basicPath5, "/", files5[1]), skip = 2), byrow = TRUE, ncol = 3),
#                dename(LoadedSimulations12[["values"]][,1], kind = "dimnames")
#   )
# })  
# 
# # differnt header form
# basicPath6 = "/home/kristina/Desktop/s4p/data/fileFormats/test_txtMulticolumn_0lineHeader"
# files6 = list.files(basicPath6)
# LoadedSimulations13 = loadSimulations(
#   basicPath = basicPath6,
#   savePath = "/home/kristina/Desktop/s4p_article/test/t13_",
#   overwrite = TRUE,
#   readBy = "scan",
#   sep = ","
# ) 
# test_that("loadSimulations txt kind multilayer 0line header",{
#   expect_equal(nrow(LoadedSimulations13[["values"]]), 100)
#   expect_equal(dename(LoadedSimulations11[["values"]][,1], kind = "dimnames"),
#                dename(LoadedSimulations13[["values"]][,1], kind = "dimnames")
#   )
# })  # 'ugly' layer names (as expected)

# number of columns varies
basicPath4a = "/home/kristina/Desktop/s4p/data/fileFormats/text_multicolumn_plumes"
expect_error(
  LoadedSimulations11a <- loadSimulations(
  basicPath = basicPath4a,
  savePath = "/home/kristina/Desktop/s4p_article/test/t11a_",
  overwrite = TRUE,
  readBy = "scan",
  skip = 1,
  sep = ",")
)  
expect_error(
  LoadedSimulations11b <- loadSimulations(
    basicPath = basicPath4a,
    savePath = "/home/kristina/Desktop/s4p_article/test/t11b_",
    overwrite = TRUE,
    readBy = "scan",
    skip = 1,
    sep = ",",
    region = region)
)  


# B)               c)
LoadedSimulations14 = loadSimulations(# interpret kinds as plumes
  basicPath = basicPath4,
  savePath = "/home/kristina/Desktop/s4p_article/test/t14_",
  multilayer = "plumes",
  overwrite = TRUE,
  readBy = "scan",
  skip = 1
)               

test_that("loadSimulations txt several plumes in files",{
  # correct dimensions
  expect_equal(nrow(LoadedSimulations14[["values"]]), 100)
  expect_equal(ncol(LoadedSimulations14[["values"]]), 15)
  expect_equal(nlayers(LoadedSimulations14[["values"]]), 1)     
  # correct values (layer order as in loadedSimulations2)
  i = sample.int(5,1); j = sample.int(3,1)
  expect_equal(LoadedSimulations14[["values"]][,(i-1) * 3 + j][,1],
               LoadedSimulations11[["values"]][,i][,j])
})

# number of columns differs among files
LoadedSimulations14a <- loadSimulations(
  basicPath = basicPath4a,
  savePath = "/home/kristina/Desktop/s4p_article/test/t14a_",
  multilayer = "plumes",
  overwrite = TRUE,
  readBy = "scan",
  skip = 1, sep = ",")
test_that("loadSimulations txt several plumes in files of differing layer number",{
  expect_equal(ncol(LoadedSimulations14a[["values"]]), 5)
})

  

# B)               a)                           i)
LoadedSimulations15 = loadSimulations(
  basicPath = basicPath3,
  savePath = "/home/kristina/Desktop/s4p_article/test/t15_",
  overwrite = TRUE,
  region = region
)

test_that("loadSimulations txt kind in files, region",{
  # correct dimensions
  expect_equal(nrow(LoadedSimulations15[["values"]]), sum(region))
  expect_equal(ncol(LoadedSimulations15[["values"]]), 5)
  expect_equal(nlayers(LoadedSimulations15[["values"]]), 3)    
  # correct values
  i = sample.int(5,1); j = sample.int(3,1)
  expect_equal(LoadedSimulations15[["values"]][,i][,j],
               LoadedSimulations10[["values"]][,i][,j][region]
               )
})

test_that("loadSimulations txt several plumes in files of differing layer number",{
  expect_equal(ncol(LoadedSimulations14a[["values"]]), 5)
})

# B)               b)                           i)
LoadedSimulations16 = loadSimulations(
  basicPath = basicPath4,
  savePath = "/home/kristina/Desktop/s4p_article/test/t16_",
  overwrite = TRUE,
  readBy = "scan",
  skip = 1,
  region = region
) 

test_that("loadSimulations txt kind multilayer, region",{
  # correct dimensions
  expect_equal(nrow(LoadedSimulations16[["values"]]), sum(region))

  expect_equal(dename(getValues(LoadedSimulations16[["values"]]), kind = "dimnames"),
               dename(getValues(LoadedSimulations15[["values"]]), kind = "dimnames")
  )
})  
# number of columns differs among files
expect_error(
  LoadedSimulations16a <- loadSimulations(
    basicPath = basicPath4a,
    savePath = "/home/kristina/Desktop/s4p_article/test/t16a_",
    region = region,
    overwrite = TRUE,
    readBy = "scan",
    skip = 1, sep = ",")
  )

# B)               c)                           i)
LoadedSimulations17 = loadSimulations(# interpret kinds as plumes
  basicPath = basicPath4,
  savePath = "/home/kristina/Desktop/s4p_article/test/t17_",
  multilayer = "plumes",
  overwrite = TRUE,
  readBy = "scan",
  skip = 1,
  region = region
)               
test_that("loadSimulations txt several plumes in files, region",{
  # correct dimensions
  expect_equal(nrow(LoadedSimulations17[["values"]]), sum(region))
  expect_equal(ncol(LoadedSimulations17[["values"]]), 15)     
  # correct values
  i = sample.int(5,1); j = sample.int(3,1)
  expect_equal(as.numeric(LoadedSimulations17[["values"]][,(i-1) * 3 + j]),
               LoadedSimulations16[["values"]][[j]][,i])
})

# number of columns differs among files
LoadedSimulations17a <- loadSimulations(
    basicPath = basicPath4a,
    savePath = "/home/kristina/Desktop/s4p_article/test/t17a_",
    region = region,
    overwrite = TRUE,
    multilayer = "plumes",
    readBy = "scan",
    skip = 1, sep = ",")
test_that("loadSimulations txt several plumes in files of differing layer number, region",{
  expect_equal(ncol(LoadedSimulations17a[["values"]]), 5)
  expect_equal(nrow(LoadedSimulations17a[["values"]]), sum(region))
})


#-------------
if(FALSE){
# examples with big data (= generate example data)
time_loadSimulationsLocal = system.time(
  RadioactivePlumes_local <- loadSimulations(
    basicPath = "/home/kristina/Desktop/s4p_article/data_shortrange/data_shortrange",
    savePath = "/home/kristina/Desktop/s4p/data/radioactivePlumesLocal_",
    overwrite = TRUE,
    readBy = "scan",
    skip = 2)  
  )
#    user  system elapsed 
#  42.862   1.907  46.063
 save(RadioactivePlumes_local,
      file = "/home/kristina/Desktop/s4p/data/radioactivePlumesLocal.rda")

time_loadSimulationsArea = system.time(
  RadioactivePlumes_area <- loadSimulations(
    basicPath = "/home/kristina/Desktop/s4p_article/data_longrange",
    savePath = "/home/kristina/Desktop/s4p/data/radioactivePlumesArea_",
    overwrite = TRUE,
    bBox = c(3476000, 3976000, 2480000, 2680000))  
) 
#      user   system  elapsed 
#  7598.298   90.096 8711.413 # this number refers to a run with a previous version
save(RadioactivePlumes_area,
     file = "/home/kristina/Desktop/s4p/data/radioactivePlumesArea.rda")
}

