##################################################################
# test simulationsApply                                          #
##################################################################

library(testthat)
test_that("simulationsApply",
{
  library(s4p)
  data(radioactivePlumes_local)
  source("/home/kristina/Desktop/s4p/R/simulationsApply.R")
  source("/home/kristina/Desktop/s4p/R/replaceDefault.R")
  radioactivePlumes_local@plumes$totalDose = summaryPlumes(radioactivePlumes_local, fun = sum)[["summaryPlumes"]]
  radioactivePlumes_local@plumes$maxDose = summaryPlumes(radioactivePlumes_local, fun = max, values = "maxdose")[["summaryPlumes"]]
  
  i = sample(nLocations(radioactivePlumes_local),1)
  j = sample(nPlumes(radioactivePlumes_local),1)
  
  ### fun
  #### fun fails check (no function)
  fun0 = 1
  simulationsApply_fun0 <- simulationsApply( # no warning (warning would be nice)
    simulations = radioactivePlumes_local,
    fun = fun0,
    chunksize = 3e+7 # needed, else chunks
  )  
  expect_equal(
    simulationsApply_fun0,
    list()
  )
  rm(simulationsApply_fun0)
  #### nout == 1
   fun1 = function(x, nout = 1){
     max(x, na.rm = TRUE)
   }
  #### no chunks
   simulationsApply_fun1 = simulationsApply(
     simulations = radioactivePlumes_local,
     fun = fun1,
     chunksize = 3e+7, # needed, else chunks
     keepSubset = TRUE
   )
  expect_equivalent(simulationsApply_fun1[["subset"]][i,,], 
                    radioactivePlumes_local@values[i,]
  )
  expect_equivalent(simulationsApply_fun1[["subset"]][,j,], 
                    radioactivePlumes_local@values[,j]
  )
  expect_equivalent(simulationsApply_fun1[["result_global"]], 
                    max(radioactivePlumes_local@values[,], na.rm = TRUE))
  #### without subset
  simulationsApply_fun2 = simulationsApply(
    simulations = radioactivePlumes_local,
    fun = fun1,
    chunksize = 3e+7 
  )
  expect_equivalent(simulationsApply_fun1[[1]], 
                    simulationsApply_fun2[[1]]      
  )
  expect_equal(
    length(simulationsApply_fun2), 1)
  
  #### chunks  -> fun cannot be applied as not all data loaded at once -> result empty
  expect_warning(
    simulationsApply_fun2a <- simulationsApply(
      simulations = radioactivePlumes_local,
      fun = fun1, chunksize = 1e+6
    )  
  )
  expect_equal(simulationsApply_fun2a, 
               list()
  )
  rm(simulationsApply_fun2a)   

  #### nout > 1
   fun2 = function(x, nout = 2){
     y = x[,,1]
     z = c(max(y, na.rm = TRUE), min(y, na.rm = TRUE))
     return(z)
   }
   simulationsApply_fun3 = simulationsApply(
     simulations = radioactivePlumes_local,
     fun = fun2, 
     chunksize = 3e+7
   ) 
  radioactivePlumes_local1 = 
    subset(radioactivePlumes_local@values, 1)[,]# as layers are the third dimension, this is how to select them
  expect_equivalent(
    simulationsApply_fun3[["result_global"]],
    c(max(radioactivePlumes_local1, na.rm = TRUE),
      min(radioactivePlumes_local1, na.rm = TRUE))
  )  
  
  # wrong number of nout: does not matter here; should actually give an error
  fun3 = function(x, nout = 2){
    max(x, na.rm = TRUE)
  }
  simulationsApply_fun4 = simulationsApply(
    simulations = radioactivePlumes_local,
    fun = fun3, 
    chunksize = 3e+7
  )
  expect_equal(
    simulationsApply_fun2, 
    simulationsApply_fun4
  )
  rm(simulationsApply_fun2)   
  rm(simulationsApply_fun4)

  ### fun_p
  #### fun fails check (extra parameter without default)
  fun_p3a = function(x, nout = 1, weight){
    sum(x[,1] * weight, na.rm = TRUE)
  }
  expect_warning(
    simulationsApply_p1a <- simulationsApply(
      simulations = radioactivePlumes_local,
      fun_p = fun_p3a
    )    
  )
  expect_equal(
    simulationsApply_p1a,
    list()
  )
  rm(simulationsApply_p1a)
  
  #### nout = 1
  fun_p3 = function(x, nout = 1){
    sum(x[,1], na.rm = TRUE)
  }
  simulationsApply_p1 = simulationsApply(
    simulations = radioactivePlumes_local,
    fun_p = fun_p3
  )
  simulationsApply_p2 = simulationsApply(
    simulations = radioactivePlumes_local,
    fun_p = fun_p3,
    chunksize = 5e+6
  )
  expect_equal(
    simulationsApply_p2,
    simulationsApply_p1
  ) 
  rm(simulationsApply_p2) 

  simulationsApply_p3 = simulationsApply(
    simulations = radioactivePlumes_local,
    fun_p = fun_p3,
    keepSubset = TRUE
  )
  expect_equal(
    simulationsApply_p3[[1]],
    simulationsApply_p1[[1]]
  )
  #### subset
  expect_equal(
    names(simulationsApply_p3),
    c(names(simulationsApply_p1), "subset")
  )
  rm(simulationsApply_p3)

  #### no subset, as too big (according to chunksize)
  expect_warning(
    simulationsApply_p4 <- simulationsApply(
      simulations = radioactivePlumes_local,
      fun_p = fun_p3,
      chunksize = 5e+6 
    )  
  )
  expect_equal(
    simulationsApply_p4,
    simulationsApply_p1
  )
  rm(simulationsApply_p1)
  rm(simulationsApply_p4)
  
  #### nout > 1
  fun_p1 = function(x, nout = 3){
    x_df = data.frame(maxdose = x[,2], finaldose = x[,1])
    lm_x = lm(finaldose ~ maxdose, x_df)
    out = c(lm_x[[1]], Rsq = summary(lm_x)[[8]])
  }
  
  simulationsApply_p5 = simulationsApply(
    simulations = radioactivePlumes_local,
    fun_p = fun_p1,
    keepSubset = TRUE
  )
  plot(simulationsApply_p5[["subset"]][,j,2], simulationsApply_p5[["subset"]][,j,1])
  abline(a = simulationsApply_p5[["result_plumes"]][j,"(Intercept)"], 
         b = simulationsApply_p5[["result_plumes"]][j,"maxdose"])

  lm_j = lm(finaldose ~ maxdose, as.data.frame(radioactivePlumes_local@values[,j]))
  expect_equivalent(simulationsApply_p5[["subset"]],
                    simulationsApply_fun1[["subset"]]
  ) 
  expect_equivalent(simulationsApply_p5[["result_plumes"]][j,],
                    c(lm_j[[1]], summary(lm_j)[[8]])
  ) 
  simulationsApply_p6 = simulationsApply(
    simulations = radioactivePlumes_local,
    fun_p = fun_p1,
    chunksize = 5e+6
  )
  expect_equivalent(simulationsApply_p5[["result_plumes"]],
                    simulationsApply_p6[["result_plumes"]])
  rm(simulationsApply_p5) 
  rm(simulationsApply_p6)

  fun_p2 = function(x, nout = 4){# wrong nout
    x_df = data.frame(maxdose = x[,2], finaldose = x[,1])
    lm_x = lm(finaldose ~ maxdose, x_df)
    out = c(lm_x[[1]], Rsq = summary(lm_x)[[8]])
  }
  expect_warning(
    simulationsApply_p6a <- simulationsApply(# wrong nout generates warning (and empty output) when processed in chunks
    simulations = radioactivePlumes_local,
    fun_p = fun_p2,
    chunksize = 5e+6)
  )  
  expect_equivalent(
    simulationsApply_p6a,
    list()
  )
  rm(simulationsApply_p6a)
  
  ### fun_p, fun_Rp
  fun_Rp1 = function(x, nout = 1, weight){
    max(x * weight[,"totalDose"], na.rm = TRUE)
  }
  expect_warning( # use of fun_Rp without fun_p
    simulationsApply_p6b <- simulationsApply(
      simulations = radioactivePlumes_local,
      fun_Rp = fun_Rp1
    )
  )
  expect_equivalent(
    simulationsApply_p6b,
    list()
  )
  rm(simulationsApply_p6b)   

#   #### nout == 1  
   detection = function(x, threshold = 1e-7, nout = 1){ # for all plumes determine how many sensors would detect them
     sum(x[,2] > threshold, na.rm = TRUE)
   }
  fun_Rp3a = function(x, nout = 1){# fails check (no 'weight') -> not applied
    mean((x >= 2) * weight[,"totalDose"], na.rm = TRUE)
  }
  expect_warning(
    simulationsApply_pRp1a <- simulationsApply(
      simulations = radioactivePlumes_local,
      fun_p = detection,
      fun_Rp = fun_Rp3a,
      keepSubset = TRUE
    )    
  )
  rm(simulationsApply_pRp1)   

   fun_Rp3 = function(x, nout = 1, weight){# how many of the plumes were detected by at least 2 sensors?
     mean((x >= 2) * weight[,"totalDose"], na.rm = TRUE)
   }
   simulationsApply_pRp1 = simulationsApply(
     simulations = radioactivePlumes_local,
     fun_p = detection,
     fun_Rp = fun_Rp3,
     keepSubset = TRUE
   )
  expect_equal(
   simulationsApply_pRp1[c(1,3)],
   simulationsApply_pRp1a # result of fun_Rp missing
  )
  rm(simulationsApply_pRp1a)
  expect_equivalent(simulationsApply_pRp1[["result_plumes"]][j],
                    sum(simulationsApply_pRp1[["subset"]][,j,2] > 1e-7, na.rm = TRUE)
  )
  simulationsApply_pRp2 = simulationsApply(
    simulations = radioactivePlumes_local,
    fun_p = detection,
    fun_Rp = fun_Rp3,
    chunksize = 5e+6
  )  
  expect_equivalent(simulationsApply_pRp1[["result_global_plumes"]],
                    simulationsApply_pRp2[["result_global_plumes"]]) 
  expect_equivalent(simulationsApply_pRp1[["result_plumes"]],     
                    simulationsApply_pRp2[["result_plumes"]]) 
  rm(simulationsApply_pRp2)
  
#### nout > 1
   fun_p2 = function(x, nout = 3){# for all plumes fit a linear model finaldose ~ maxdose
     x_df = data.frame(maxdose = x[,2], finaldose = x[,1])
     lm_x = lm(finaldose ~ maxdose, x_df)
     out = c(lm_x[[1]], Rsq = summary(lm_x)[[8]])
   }
   fun_Rp2 = function(x, nout = 2, weight){# how many of the models have Rsq < 0.5 (and slope <1)?
     c(mean(x[,3] <= 0.5, na.rm = TRUE),
       mean((x[,2] <= 1) * weight[,"totalDose"], na.rm = TRUE))
   }
   simulationsApply_pRp3 = simulationsApply(
     simulations = radioactivePlumes_local,
     fun_p = fun_p2,
     fun_Rp = fun_Rp2
   )
  expect_equivalent(simulationsApply_pRp3[["result_global_plumes"]][1],
                    mean(simulationsApply_pRp3[["result_plumes"]][,"Rsq"] <= 0.5, na.rm = TRUE)
  )

  fun_p2a = function(nout = 3){# no 'x'-> fails?
    x_df = data.frame(maxdose = x[,2], finaldose = x[,1])
    lm_x = lm(finaldose ~ maxdose, x_df)
    out = c(lm_x[[1]], Rsq = summary(lm_x)[[8]])
  }
  expect_warning(# fun_p fails -> no result of fun_Rp
    simulationsApply_pRp3a <- simulationsApply(
      simulations = radioactivePlumes_local,
      fun_p = fun_p2a,
      fun_Rp = fun_Rp2
    ) 
  )
  expect_equal(
    simulationsApply_pRp3a,
    list()
  )
  rm(simulationsApply_pRp3a)

 simulationsApply_pRp4 = simulationsApply(
   simulations = radioactivePlumes_local,
   fun_p = fun_p2,
   fun_Rp = fun_Rp2,
   chunksize = 5e+6
 )  
 expect_equivalent(simulationsApply_pRp3[["result_global_plumes"]],
                   simulationsApply_pRp4[["result_global_plumes"]])
 expect_equivalent(simulationsApply_pRp3[["result_plumes"]],
                   simulationsApply_pRp4[["result_plumes"]]
 )
 rm(simulationsApply_pRp4); gc()

  ### fun_l, fun_Rl
  #### fun_l fails check -> fun_Rl as well not applied
  fun_l1a = function(x, nout = 1, weight){
    excMeasured = sum(x[,2] > 1e-3, na.rm = TRUE)
  }
   fun_Rl1 = function(x, nout = 1, weight){
     mean(x * weight[,"index"], na.rm = TRUE)
   }
  expect_warning(
    simulationsApply_lRl1a <- simulationsApply(
      simulations = radioactivePlumes_local,
      fun_l = fun_l1a,
      fun_Rl = fun_Rl1,
      keepSubset = TRUE
    )    
  )
  rm(simulationsApply_lRl1)

  #### nout = 1
   fun_l1 = function(x, nout = 1){ # how many threshold exceedances are measured here?
     excMeasured = sum(x[,2] > 1e-3, na.rm = TRUE)
   }
 
   simulationsApply_lRl1 = simulationsApply(
     simulations = radioactivePlumes_local,
     fun_l = fun_l1,
     fun_Rl = fun_Rl1,
     keepSubset = TRUE
   )
  expect_equivalent(sum(simulationsApply_lRl1[["subset"]][i,,2] > 1e-3), 
                    simulationsApply_lRl1[["result_locations"]][i]
  )
  expect_equivalent(mean(simulationsApply_lRl1[["result_locations"]] * 
                           radioactivePlumes_local@locations@data[,"index"], na.rm = TRUE),
                    simulationsApply_lRl1[["result_global_locations"]]
  )
  expect_equal(
    simulationsApply_lRl1["subset"],
    simulationsApply_lRl1a
  )
  rm(simulationsApply_lRl1a)

  simulationsApply_lRl2 = simulationsApply(
    simulations = radioactivePlumes_local,
    fun_l = fun_l1,
    fun_Rl = fun_Rl1,
    chunksize = 5e+6
  )
  expect_equivalent(simulationsApply_lRl1[[1:2]],  
                  simulationsApply_lRl2[[1:2]]
  )
  rm(simulationsApply_lRl2)

  #### nout > 1
   fun_l2 = function(x, nout = 2){ # how many threshold exceedances are measured here?
     excMeasured = c(sum(x[,2] > 1e-3, na.rm = TRUE),
                     sum(x[,2] > 1e-5, na.rm = TRUE))
   }
   fun_Rl2 = function(x, nout = 2, weight, y = 2 * x){# extra parameter, with default: works
     c(0.5 * mean(y[,1], na.rm = TRUE),
       mean(x[,1]/weight[,"index"], na.rm = TRUE))
   }
   simulationsApply_lRl3 = simulationsApply(
     simulations = radioactivePlumes_local,
     fun_l = fun_l2,
     fun_Rl = fun_Rl2
   )
  simulationsApply_lRl4 = simulationsApply(
    simulations = radioactivePlumes_local,
    fun_l = fun_l2,
    fun_Rl = fun_Rl2,
    chunksize = 5e+6
  )
  expect_equivalent(mean(simulationsApply_lRl3[["result_locations"]][,1]/radioactivePlumes_local@locations@data[,"index"], na.rm = TRUE), 
                    simulationsApply_lRl3[["result_global_locations"]][2]
  )
  expect_equivalent(simulationsApply_lRl3[["result_global_locations"]],
                    simulationsApply_lRl4[["result_global_locations"]])
  expect_equivalent(simulationsApply_lRl3[["result_locations"]],
                    simulationsApply_lRl4[["result_locations"]])  
  fun_Rl2a = function(x, nout = 2, weight, ...){# extra parameter -> fails
    c(mean(x[,1], na.rm = TRUE),
    mean(x[,2]/weight[,"index"], na.rm = TRUE))
  }
  expect_warning(
    simulationsApply_lRl3a <- simulationsApply(
      simulations = radioactivePlumes_local,
      fun_l = fun_l2,
      fun_Rl = fun_Rl2a # fails -> result missing
    )        
  )  
  expect_equivalent(simulationsApply_lRl3[["result_locations"]],
                    simulationsApply_lRl3a[["result_locations"]]) 
  
  ### fun_pl, fun_R_pl
  #### errors in fun
  detection2a = function(x, threshold = 1e-7, nout = 1, ...){ 
    sum(x[2] > threshold, na.rm = TRUE)
  }
   fun_Rpl1 = function(x, nout = 1, weight_l = 10, weight_p = 10){
     mean(x * weight_l[,"index"], na.rm = TRUE)
   }
  expect_warning( # wrong param in detection2a -> no results
    simulationsApply_plRpl1a <- simulationsApply( 
      radioactivePlumes_local,
      fun_pl = detection2a,
      fun_Rpl = fun_Rpl1
    )     
  )
  expect_equal(
    simulationsApply_plRpl1a,
    list()
  )
  rm(simulationsApply_plRpl1a)

  #### errors in fun_Rpl
   detection2 = function(x, threshold = 1e-7, nout = 1){ 
     sum(x[2] > threshold, na.rm = TRUE)
   }
  fun_Rpl1a = function(x, nout = 1, weight_l = 10, weight_p = 10, ...){
    mean(x * weight_p[,"totalDose"], na.rm = TRUE)
  }
  expect_warning( # wrong param in fun_Rpl1a -> no result of it
    simulationsApply_plRpl1b <- simulationsApply( 
      radioactivePlumes_local,
      fun_pl = detection2,
      fun_Rpl = fun_Rpl1a
    )     
  )

  #### nout == 1
   simulationsApply_plRpl1 = simulationsApply( 
       radioactivePlumes_local,
       fun_pl = detection2,
       fun_Rpl = fun_Rpl1
   )
  expect_equivalent(
    simulationsApply_plRpl1b,
    simulationsApply_plRpl1[1]
  )
  rm(simulationsApply_plRpl1b)

  expect_equivalent(
    simulationsApply_plRpl1[["result_global_locationsplumes"]],
    mean(radioactivePlumes_local@locations@data[,"index"] * 
           getValues(simulationsApply_plRpl1[["result_locationsplumes"]]), na.rm = TRUE) 
  )

  #### chunks
  expect_warning(# result of fun_pl not in memory, no file -> empty output
    simulationsApply_pl1c <- simulationsApply( 
      radioactivePlumes_local,
      fun_pl = detection2,
      fun_Rpl = fun_Rpl1,
      chunksize = 5e+6,
      nameSave = FALSE
    )
  )
  expect_equal(
    simulationsApply_pl1c,
    list()
  )
  rm(simulationsApply_pl1c)

  expect_warning(# result of fun_pl not in memory -> fun_Rpl cannot be applied, files generated? 
    simulationsApply_pl1d <- simulationsApply( 
      radioactivePlumes_local,
      fun_pl = detection2,
      fun_Rpl = fun_Rpl1,
      chunksize = 5e+6,
      nameSave = "tmp_simApply1",
      overwrite = TRUE
    )
  )
  expect_equal(
    length(simulationsApply_pl1d), 1
  )
  expect_is(
    simulationsApply_pl1d[["result_locationsplumes"]],
    "RasterBrick"
  )
  rm(simulationsApply_pl1d)

   detection2b = function(x, threshold = 1e-7, nout = 2){# nout of fun_pl is wrong 
     sum(x[2] > threshold, na.rm = TRUE)
   }
  expect_warning(# nout of fun_pl wrong -> no results
    simulationsApply_pl1e <- simulationsApply( 
      radioactivePlumes_local,
      fun_pl = detection2b,
      fun_Rpl = fun_Rpl1,
      chunksize = 5e+6,
      nameSave = "tmp_simApply1e",
      overwrite = TRUE
    )
  )
  expect_equivalent(
    simulationsApply_pl1e,# named list
    list()
  )
  rm(simulationsApply_pl1e)
  expect_false(
    file.exists("tmp_simApply1e.grd")  
  )
  expect_false(
    file.exists("tmp_simApply1e.gri")  
  )

  expect_warning( # result is saved to file
    simulationsApply_pl2 <- simulationsApply( 
      radioactivePlumes_local,
      fun_pl = detection2,
      nameSave = "tmp_simApply_pl2",#"/home/kristina/Desktop/s4p_article/test/simulationsApply_17_simApply",
      chunksize = 2e+6,
      overwrite = TRUE
    )  
  )
  expect_equivalent(
    simulationsApply_pl2[["result_locationsplumes"]], 
    simulationsApply_plRpl1[["result_locationsplumes"]]
  ) 
  rm(simulationsApply_pl2) 

  #### nout > 1
   detection3 = function(x, threshold = 1e-7, nout = 2){ 
     c(sum(x[1] > threshold, na.rm = TRUE),
       sum(x[2] > threshold, na.rm = TRUE))
   }
   fun_Rpl2 = function(x, nout = 2, weight_l = 1, weight_p = 1){ 
     c(mean(x[,1] * weight_l[,"index"], na.rm = TRUE),
       mean(x[,2]/weight_p[,"totalDose"], na.rm = TRUE))
   }
   simulationsApply_plRpl3 = simulationsApply( 
     radioactivePlumes_local,
     fun_pl = detection3,
     fun_Rpl = fun_Rpl2,
     keepSubset = TRUE
   )
  expect_equivalent(
    matrix(detection3(simulationsApply_plRpl3[["subset"]][i + 20, j + 6,]), nrow = 1),
    simulationsApply_plRpl3[["result_locationsplumes"]][i + 20, j + 6]
  )
  expect_equivalent(
    simulationsApply_plRpl3[["result_global_locationsplumes"]],
    c(mean(getValues(simulationsApply_plRpl3[["result_locationsplumes"]])[,1] * 
             radioactivePlumes_local@locations@data[,"index"], na.rm = TRUE),
      mean(getValues(simulationsApply_plRpl3[["result_locationsplumes"]])[,2] / 
             radioactivePlumes_local@plumes[,"totalDose"], na.rm = TRUE))
      
    #cellStats(simulationsApply_plRpl3[["result_locationsplumes"]], "mean", na.rm = TRUE)
  )

  #### chunks
  expect_warning(# result of fun_pl not in memory -> fun_Rpl cannot be applied 
    simulationsApply_plRpl4 <- simulationsApply( 
      radioactivePlumes_local,
      fun_pl = detection3,
      fun_Rpl = fun_Rpl2,
      chunksize = 2e+6,
      nameSave = "tmp_simulationsApply_plRlp4",
      overwrite = TRUE
    )    
  )
  rm(simulationsApply_plRpl4)
  expect_warning(# result is saved to file
    simulationsApply_pl4 <- simulationsApply( 
      radioactivePlumes_local,
      fun_pl = detection3,
      chunksize = 2e+6,
      nameSave = "tmp_simulationsApply_plp4",
      overwrite = TRUE
    )   
  )
  
  expect_equivalent(
    simulationsApply_pl4[["result_locationsplumes"]],
    simulationsApply_plRpl3[["result_locationsplumes"]]
  )
  rm(simulationsApply_pl4)
  expect_true(
    file.exists("tmp_simulationsApply_plp4_locationsplumes.gri") 
  )

#   expect_equal(
#     matrix(simulationsApply_19[["result_locationsplumes"]][2,i + 1:20, j + 1:10], ncol = 10, byrow = TRUE), # ~1dim
#     simulationsApply_18[["result_locationsplumes"]][2,i + 1:20, j + 1:10]
#   ) 

# combination of p,l,pl
### chunks will cause fun & fun_Rpl to have no result
  expect_warning(
    simulationsApply_p_l_pl_R1 <- simulationsApply(
      simulations = radioactivePlumes_local,
      fun = fun1, # from sA_fun1
      fun_p = detection, # from sA_pRp1
      fun_Rp = fun_Rp3, # from sA_pRp1
      fun_l = fun_l1, # from sA_lRl1
      fun_Rl = fun_Rl1, # from sA_lRl1
      fun_pl = detection2, # from sA_plRpl1
      fun_Rpl = fun_Rpl1, # from sA_plRpl1
      chunksize = 1e+6,
      nameSave = "tmp_p_l_pl_R1",
      overwrite = TRUE
    )  
  )
  expect_equal(
    names(simulationsApply_p_l_pl_R1),
    c("result_locations", "result_global_locations","result_plumes","result_global_plumes", "result_locationsplumes")
  )
  expect_equal(
    simulationsApply_p_l_pl_R1[["result_plumes"]],
    simulationsApply_pRp1[["result_plumes"]]
  )
  expect_equal(
    simulationsApply_p_l_pl_R1[["result_global_plumes"]],
    simulationsApply_pRp1[["result_global_plumes"]]
  )
  expect_equal(
    simulationsApply_p_l_pl_R1[["result_locations"]],
    simulationsApply_lRl1[["result_locations"]]
  )
  expect_equal(
    simulationsApply_p_l_pl_R1[["result_global_locations"]],
    simulationsApply_lRl1[["result_global_locations"]]
  )
  expect_equivalent(
    simulationsApply_p_l_pl_R1[["result_locationsplumes"]],
    simulationsApply_plRpl1[["result_locationsplumes"]]
  )
  ## nout == 1
  #### no chunks
  simulationsApply_p_l_pl_R2 = simulationsApply(
    simulations = radioactivePlumes_local,
    fun = fun1, # from sA_fun1
    fun_p = detection, # from sA_pRp1
    fun_Rp = fun_Rp3, # from sA_pRp1
    fun_l = fun_l1, # from sA_lRl1
    fun_Rl = fun_Rl1, # from sA_lRl1
    fun_pl = detection2, # from sA_plRpl1
    fun_Rpl = fun_Rpl1, # from sA_plRpl1
    chunksize = 3e+7 
  )
  expect_equivalent(
    simulationsApply_p_l_pl_R2[names(simulationsApply_p_l_pl_R1)],
    simulationsApply_p_l_pl_R1
  )
  rm(simulationsApply_p_l_pl_R1); gc()
  expect_equivalent(
    simulationsApply_p_l_pl_R2[["result_global"]],
    simulationsApply_fun1[["result_global"]]
  )
  rm(simulationsApply_fun1)
  expect_equivalent(
    simulationsApply_p_l_pl_R2[["result_global_locationsplumes"]],
    simulationsApply_plRpl1[["result_global_locationsplumes"]]
  )
  rm(simulationsApply_p_l_pl_R2); gc() # else, R may break by lack of memory

  ## nout > 1 
  simulationsApply_p_l_pl_R3 = simulationsApply(
    simulations = radioactivePlumes_local,
    fun = fun2, # from sA_fun3
    fun_p = fun_p2, # from sA_pRp3
    fun_Rp = fun_Rp2, # from sA_pRp3
    fun_l = fun_l2, # from sA_lRl3
    fun_Rl = fun_Rl2, # from sA_lRl3
    fun_pl = detection3, # from sA_plRpl3
    fun_Rpl = fun_Rpl2, # from sA_plRpl3
    chunksize = 3e+7 
  )

  expect_equivalent(
    simulationsApply_p_l_pl_R3[["result_global"]],
    simulationsApply_fun3[["result_global"]]
  )
  rm(simulationsApply_fun3)
  expect_equivalent(
    simulationsApply_p_l_pl_R3[["result_plumes"]],
    simulationsApply_pRp3[["result_plumes"]]
  )
  expect_equivalent(
    simulationsApply_p_l_pl_R3[["result_global_plumes"]],
    simulationsApply_pRp3[["result_global_plumes"]]
  )
  rm(simulationsApply_pRp3)
  expect_equivalent(
    simulationsApply_p_l_pl_R3[["result_locations"]],
    simulationsApply_lRl3[["result_locations"]]
  )
  expect_equivalent(
    simulationsApply_p_l_pl_R3[["result_global_locations"]],
    simulationsApply_lRl3[["result_global_locations"]]
  )
  rm(simulationsApply_lRl3)
  expect_equivalent(
    simulationsApply_p_l_pl_R3[["result_locationsplumes"]],
    simulationsApply_plRpl3[["result_locationsplumes"]]
  )
  expect_equivalent(
    simulationsApply_p_l_pl_R3[["result_global_locationsplumes"]],
    simulationsApply_plRpl3[["result_global_locationsplumes"]]
  )
  rm(simulationsApply_plRpl3)
  rm(simulationsApply_p_l_pl_R3)

## subsetting (with multiple, switched order, invalid, NA)
## is weights are used: everything (especially subsetting), including weights 

#locations = integer(0) 
expect_warning(
  simulationsApply__1 <- simulationsApply(
    simulations = radioactivePlumes_local,
    locations = integer(0),
    fun = fun1, # from sA_fun1
    fun_p = detection, # from sA_pRp1
    fun_Rp = fun_Rp3, # from sA_pRp1
    fun_l = fun_l1, # from sA_lRl1
    fun_Rl = fun_Rl1, # from sA_lRl1
    fun_pl = detection2, # from sA_plRpl1
    fun_Rpl = fun_Rpl1, # from sA_plRpl1
    chunksize = 3e+7 
  )  
)
expect_equal(
  simulationsApply__1[["result_global"]],
  fun1(integer(0))
)
expect_equivalent(
  simulationsApply__1[["result_locations"]],
  matrix(nrow = 0, ncol = formals(fun_l1)$nout)
)

expect_equal(
  simulationsApply__1[["result_global_locations"]],
  fun_Rl1(x = simulationsApply__1[["result_locations"]], weight = radioactivePlumes_local@locations@data)
)
expect_equal(
  simulationsApply__1[["result_plumes"]],
  matrix(detection(matrix(nrow = 0, ncol = 3)), 
         nrow = nPlumes(radioactivePlumes_local), ncol = formals(detection)$nout)
)
expect_equal(
  simulationsApply__1[["result_global_plumes"]],
  fun_Rp3(simulationsApply__1[["result_plumes"]], weight = radioactivePlumes_local@plumes)
)
expect_equivalent(
  simulationsApply__1[["result_locationsplumes"]],
  array(detection2(array(dim = c(0, 1, 3))), 
         dim = c(0, nPlumes(radioactivePlumes_local), ncol = formals(detection)$nout))
)
  rm(simulationsApply__1)

# plumes = integer(0) ->
expect_warning(
  simulationsApply__2 <- simulationsApply(
    simulations = radioactivePlumes_local,
    plumes = integer(0),
    fun = fun1, # from sA_fun1
    fun_p = detection, # from sA_pRp1
    fun_Rp = fun_Rp3, # from sA_pRp1
    fun_l = fun_l1, # from sA_lRl1
    fun_Rl = fun_Rl1, # from sA_lRl1
    fun_pl = detection2, # from sA_plRpl1
    fun_Rpl = fun_Rpl1, # from sA_plRpl1
    chunksize = 3e+7 
  )  
)

expect_equal(
  simulationsApply__2[["result_global"]],
  fun1(integer(0))
)
expect_equivalent(
  simulationsApply__2[["result_locations"]],
  matrix(fun_l1(matrix(nrow = 0, ncol = 3)), 
        nrow = nLocations(radioactivePlumes_local), ncol = formals(detection)$nout)
)
expect_equal(
  simulationsApply__2[["result_global_locations"]],
  fun_Rl1(simulationsApply__2[["result_locations"]], weight = radioactivePlumes_local@locations@data)
)
expect_equivalent(
  simulationsApply__2[["result_plumes"]],
  matrix(nrow = 0, ncol = formals(detection)$nout)
)
expect_equal(
  simulationsApply__2[["result_global_plumes"]],
  fun_Rp3(simulationsApply__2[["result_plumes"]], weight = radioactivePlumes_local@plumes)
)
expect_equivalent(
  simulationsApply__2[["result_locationsplumes"]],
  array(detection2(array(dim = c(1, 0, 3))), 
        dim = c(nLocations(radioactivePlumes_local), 0, ncol = formals(detection)$nout))
)
  rm(simulationsApply__2)

# test with functions that have special output for empty input?

# test subsetting (not empty, including NA, invalid (and therefore empty))
# big, including reversed order
locationsI = sample.int(nLocations(radioactivePlumes_local), 4000)
locationsI = c(locationsI, rev(locationsI))
plumesI = sample.int(nPlumes(radioactivePlumes_local), 400)
plumesI = c(plumesI, rev(plumesI))
i = sample.int(8000, 1)
j = sample.int(800, 1)
### nout == 1, no chunks
simulationsApply__3 = simulationsApply(
  simulations = radioactivePlumes_local,
  locations = locationsI,
  plumes = plumesI,
  fun = fun1, # from sA_fun1
  fun_p = detection, # from sA_pRp1
  fun_Rp = fun_Rp3, # from sA_pRp1
  fun_l = fun_l1, # from sA_lRl1
  fun_Rl = fun_Rl1, # from sA_lRl1
  fun_pl = detection2, # from sA_plRpl1
  fun_Rpl = fun_Rpl1, # from sA_plRpl1
  chunksize = 3e+7,
  keepSubset = TRUE
)
#!#! Error in subset.Simulations: 
# invalid class “Simulations” oject: The number of 'locations' ( 4000 )
#differs from the number of rows of 'values' ( 8000 ).
#They have to agree, the rows of 'values' must belong to the locations that are represented by the rows of 'locations'.

# problem: subsetSDF deletes duplicates
# not for SPointsDF (and all other classes derived from sp, but for SIndexDF and SPolygridDF)
subset_radioactivePlumes_local =  subset(x = radioactivePlumes_local, locations = locationsI, plumes = plumesI, valuesOnly = TRUE, overwrite = TRUE)
# subset_radioactivePlumes_local =  radioactivePlumes_local@values[locationsI, plumesI]
#Subset_radioactivePlumes_local = array(dim = c(length(locationsI), length(plumesI), 3))
#for (i in 1:3){
#  Subset_radioactivePlumes_local[,,i] = matrix(subset_radioactivePlumes_local[,i], 
#                                               byrow = TRUE, nrow = length(locationsI))  
#}
#subset(radioactivePlumes_local, locations = locationsI, plumes = plumesI, saveName = "subset_simulationsApply1", overwrite = TRUE)
#expect_equal(
#  subset_radioactivePlumes_local,
#  simulationsApply__3[["subset"]]
#)
expect_equal(
  simulationsApply__3[["result_global"]],
  fun1(getValues(subset_radioactivePlumes_local))
#  fun1(simulationsApply__3[["subset"]])
)
expect_equal(
  simulationsApply__3[["result_locations"]][i],
  fun_l1(subset_radioactivePlumes_local[i,,])
#  fun_l1(simulationsApply__3[["subset"]][i,,])
)
expect_equal(
  simulationsApply__3[["result_global_locations"]],
  fun_Rl1(simulationsApply__3[["result_locations"]], 
          weight = radioactivePlumes_local@locations@data[locationsI,,drop = FALSE])
)
expect_equal(
  simulationsApply__3[["result_plumes"]][j],
  detection(subset_radioactivePlumes_local[,j,])
)
expect_equal(
  simulationsApply__3[["result_global_plumes"]],
  fun_Rp3(simulationsApply__3[["result_plumes"]],
          weight = radioactivePlumes_local@plumes[plumesI,,drop = FALSE])
)
expect_equivalent(
  as.numeric(simulationsApply__3[["result_locationsplumes"]][i,j]),
  detection2(subset_radioactivePlumes_local[i,j,])
)
expect_equivalent(
  simulationsApply__3[["result_global_locationsplumes"]],
  fun_Rpl1(getValues(simulationsApply__3[["result_locationsplumes"]]),
           weight_p = radioactivePlumes_local@plumes[plumesI,,drop = FALSE],
           weight_l = radioactivePlumes_local@locations@data[locationsI,,drop = FALSE])
)


# nout > 1, no chunks
simulationsApply__4 = simulationsApply( # nout > 1
  simulations = radioactivePlumes_local,
  locations = locationsI,
  plumes = plumesI,
  fun = fun2, # from sA_fun3
  fun_p = fun_p2, # from sA_pRp3
  fun_Rp = fun_Rp2, # from sA_pRp3
  fun_l = fun_l2, # from sA_lRl3
  fun_Rl = fun_Rl2, # from sA_lRl3
  fun_pl = detection3, # from sA_plRpl3
  fun_Rpl = fun_Rpl2, # from sA_plRpl3
  chunksize = 3e+7,
  keepSubset = TRUE
)
expect_equal(
  simulationsApply__4[["result_global"]],
  fun2(simulationsApply__4[["subset"]])
)
expect_equal(
  simulationsApply__4[["result_locations"]][i,],
  fun_l2(simulationsApply__4[["subset"]][i,,])
)
expect_equal(
  simulationsApply__4[["result_global_locations"]],
  fun_Rl2(simulationsApply__4[["result_locations"]], weight = radioactivePlumes_local@locations@data[locationsI,,drop = FALSE])
)
expect_equal(
  simulationsApply__4[["result_plumes"]][j,],
  fun_p2(simulationsApply__4[["subset"]][,j,])
)
expect_equal(
  simulationsApply__4[["result_global_plumes"]],
  fun_Rp2(simulationsApply__4[["result_plumes"]], weight = radioactivePlumes_local@plumes[plumesI,,drop = FALSE])
)
expect_equivalent(
  as.numeric(simulationsApply__4[["result_locationsplumes"]][i,j,]),
  detection3(simulationsApply__4[["subset"]][i,j,])
)
expect_equivalent(
  simulationsApply__4[["result_global_locationsplumes"]],
  fun_Rpl2(getValues(simulationsApply__4[["result_locationsplumes"]]),
           weight_p = radioactivePlumes_local@plumes[plumesI,,drop = FALSE],
           weight_l = radioactivePlumes_local@locations@data[locationsI,,drop = FALSE])
)

# nout == 1,  chunks
simulationsApply__5 = simulationsApply(
  simulations = radioactivePlumes_local,
  locations = locationsI,
  plumes = plumesI,
  fun = fun1, # from sA_fun1
  fun_p = detection, # from sA_pRp1
  fun_Rp = fun_Rp3, # from sA_pRp1
  fun_l = fun_l1, # from sA_lRl1
  fun_Rl = fun_Rl1, # from sA_lRl1
  fun_pl = detection2, # from sA_plRpl1
  fun_Rpl = fun_Rpl1, # from sA_plRpl1
  chunksize = 5e+6,
  nameSave = "test_sA5",
  overwrite = TRUE
)
expect_equal(
  simulationsApply__5[1:4],
  simulationsApply__3[c(2,3,4,5)]
)
expect_equivalent(
  simulationsApply__5[[5]],
  simulationsApply__3[[6]]
)
rm(simulationsApply__3)


# nout > 1, chunks
simulationsApply__6 = simulationsApply( # nout > 1
  simulations = radioactivePlumes_local,
  locations = locationsI,
  plumes = plumesI,
  fun = fun2, # from sA_fun3
  fun_p = fun_p2, # from sA_pRp3
  fun_Rp = fun_Rp2, # from sA_pRp3
  fun_l = fun_l2, # from sA_lRl3
  fun_Rl = fun_Rl2, # from sA_lRl3
  fun_pl = detection3, # from sA_plRpl3
  fun_Rpl = fun_Rpl2, # from sA_plRpl3
  chunksize = 5e+6,
  overwrite = TRUE
)
expect_equal(
  simulationsApply__6[[1]],
  simulationsApply__4[[2]]
)
expect_equal(
  simulationsApply__6[[2]],
  simulationsApply__4[[3]]
)
expect_equivalent(
  simulationsApply__6[[3]],
  simulationsApply__4[[4]]
)
expect_equal(
  simulationsApply__6[[4]],
  simulationsApply__4[[5]]
)
expect_equivalent(
  simulationsApply__6[[5]],
  simulationsApply__4[[6]]
)
rm(simulationsApply__4)

# test weight
## generate values to use as weight
radioactivePlumes_local@plumes$totalDose = 
  summaryPlumes(radioactivePlumes_local, fun = sum, kinds = 1)[[2]]
radioactivePlumes_local@plumes$maxDose = 
  summaryPlumes(radioactivePlumes_local, fun = max, kinds = 2)[[2]]
radioactivePlumes_local@locations@data$totaldose = 
  summaryLocations(radioactivePlumes_local, fun = sum, kinds = 1)[[2]]
radioactivePlumes_local@locations@data$maxdose = 
  summaryLocations(radioactivePlumes_local, fun = max, kinds = 2)[[2]]

## functions with nout = 1, nout > 1 that use weights
## nout == 1
fun_Rp4 = function(x, nout = 1, weight){
  mean(x >= 2 * weight$totalDose, na.rm = TRUE)
}
fun_Rl3 = function(x, nout = 1, weight){
  mean(x * weight$maxdose, na.rm = TRUE)
}
simulationsApply__7 = simulationsApply(# nout = 1, chunks
  simulations = radioactivePlumes_local,
  locations = locationsI,
  plumes = plumesI,
  fun = fun1, # from sA_fun1
  fun_p = detection, # from sA_pRp1
  fun_Rp = fun_Rp4, # from sA_pRp1
  fun_l = fun_l1, # from sA_lRl1
  fun_Rl = fun_Rl3, # from sA_lRl1
  fun_pl = detection2, # from sA_plRpl1
  fun_Rpl = fun_Rpl1, # from sA_plRpl1
  chunksize = 5e+6,
  nameSave = "test_sA7",
  overwrite = TRUE
)
expect_equivalent(
  simulationsApply__7[c(1,3,5)],
  simulationsApply__5[c(1,3,5)]
)
rm(simulationsApply__5)
expect_equal(
  simulationsApply__7[[2]],
  fun_Rl3(simulationsApply__7[[1]], weight = radioactivePlumes_local@locations@data[locationsI,])
)
expect_equal(
  simulationsApply__7[[4]],
  fun_Rp4(simulationsApply__7[[3]], weight = radioactivePlumes_local@plumes[plumesI,])
)

simulationsApply__7a = simulationsApply(# nout = 1, no chunks
  simulations = radioactivePlumes_local,
  locations = locationsI,
  plumes = plumesI,
  fun = fun1, # from sA_fun1
  fun_p = detection, # from sA_pRp1
  fun_Rp = fun_Rp4, # from sA_pRp1
  fun_l = fun_l1, # from sA_lRl1
  fun_Rl = fun_Rl3, # from sA_lRl1
  fun_pl = detection2, # from sA_plRpl1
  fun_Rpl = fun_Rpl1, # from sA_plRpl1
  chunksize = 3e+7
)
expect_equivalent(
  simulationsApply__7,
  simulationsApply__7a[2:6]
)
rm(simulationsApply__7, simulationsApply__7a); gc()

## nout >1
fun_Rp5 = function(x, nout = 2, weight){# how many of the models have Rsq < 0.5 (and slope <1)?
  c(mean(x[,3] <= 0.5 * weight$maxDose, na.rm = TRUE),
    mean(x[,2] <= 1 * weight$totalDose, na.rm = TRUE))
}
fun_Rl4 = function(x, nout = 2, weight, y = 2 * x){# extra parameter, with default: works
  c(0.5 * mean(y[,1] * weight$maxdose, na.rm = TRUE),
    mean(x[,1] * weight$totaldose , na.rm = TRUE))
}
simulationsApply__8 = simulationsApply( # nout > 1, chunks
  simulations = radioactivePlumes_local,
  locations = locationsI,
  plumes = plumesI,
  fun = fun2, # from sA_fun3
  fun_p = fun_p2, # from sA_pRp3
  fun_Rp = fun_Rp5, # from sA_pRp3
  fun_l = fun_l2, # from sA_lRl3
  fun_Rl = fun_Rl2, # from sA_lRl3
  fun_pl = detection3, # from sA_plRpl3
  fun_Rpl = fun_Rpl2, # from sA_plRpl3
  chunksize = 5e+6,
  nameSave = "test_sA8",
  overwrite = TRUE
)
expect_equivalent(
  simulationsApply__8[c(1,3,5)],
  simulationsApply__6[c(1,3,5)]
)
rm(simulationsApply__6)
expect_equal(
  simulationsApply__8[[2]],
  fun_Rl2(simulationsApply__8[[1]], weight = radioactivePlumes_local@locations@data[locationsI,])
)
expect_equal(
  simulationsApply__8[[4]],
  fun_Rp5(simulationsApply__8[[3]], weight = radioactivePlumes_local@plumes[plumesI,])
)

simulationsApply__8a = simulationsApply( # nout > 1, no chunks
  simulations = radioactivePlumes_local,
  locations = locationsI,
  plumes = plumesI,
  fun = fun2, # from sA_fun3
  fun_p = fun_p2, # from sA_pRp3
  fun_Rp = fun_Rp5, # from sA_pRp3
  fun_l = fun_l2, # from sA_lRl3
  fun_Rl = fun_Rl2, # from sA_lRl3
  fun_pl = detection3, # from sA_plRpl3
  fun_Rpl = fun_Rpl2, # from sA_plRpl3
  chunksize = 3e+7
)
expect_equivalent(
  simulationsApply__8,
  simulationsApply__8a[2:6]
)
rm(simulationsApply__8)

# cellStats
expect_warning(
  simulationsApply__9 <- simulationsApply(
    simulations = radioactivePlumes_local,
    locations = locationsI,
    plumes = plumesI,
    fun = fun1, # from sA_fun1
    fun_p = detection, # from sA_pRp1
    fun_Rp = fun_Rp4, # from sA_pRp1
    fun_l = fun_l1, # from sA_lRl1
    fun_Rl = fun_Rl3, # from sA_lRl1
    fun_pl = detection2, # from sA_plRpl1
    fun_Rpl = fun_Rpl1, # from sA_plRpl1
    fun_Rpl_cellStats = "mean",
    chunksize = 5e+6,
    nameSave = "test_sA9",
    overwrite = TRUE
  )  
)
expect_equal(
  simulationsApply__9[["cellStats_global_locationsplumes"]],
  cellStats(x = simulationsApply__9[["result_locationsplumes"]], "mean", asSample = FALSE)
)
rm(simulationsApply__9)
expect_warning(
  simulationsApply__10 <- simulationsApply( # nout > 1, chunks
    simulations = radioactivePlumes_local,
    locations = locationsI,
    plumes = plumesI,
    fun = fun2, # from sA_fun3
    fun_p = fun_p2, # from sA_pRp3
    fun_Rp = fun_Rp5, # from sA_pRp3
    fun_l = fun_l2, # from sA_lRl3
    fun_Rl = fun_Rl2, # from sA_lRl3
    fun_pl = detection3, # from sA_plRpl3
    fun_Rpl = fun_Rpl2, # from sA_plRpl3
    fun_Rpl_cellStats = "mean",
    chunksize = 5e+6,
    nameSave = "test_sA10",
    overwrite = TRUE
  )
)
expect_equal(
  simulationsApply__10[["cellStats_global_locationsplumes"]],
  cellStats(x = simulationsApply__10[["result_locationsplumes"]], "mean", asSample = FALSE)
)
rm(simulationsApply__10)

# simulations is raster: fun_Rp, fun_Rl need parameter 'weight', must not use it unless it has default value
## weight not used
fun_Rp6 = function(x, nout = 2, weight){# how many of the models have Rsq < 0.5 (and slope <1)?
  c(mean(x[,3] <= 0.5, na.rm = TRUE),
    mean(x[,2] <= 1, na.rm = TRUE))
}
fun_Rl5 = function(x, nout = 2, y = 2 * x, weight){# extra parameter, with default: works
  c(0.5 * mean(y[,1], na.rm = TRUE),
    mean(x[,1], na.rm = TRUE))
}
simulationsApply__11 <- simulationsApply( # nout > 1, chunks
  simulations = radioactivePlumes_local,
  locations = locationsI,
  plumes = plumesI,
  fun = fun2, # from sA_fun3
  fun_p = fun_p2, # from sA_pRp3
  fun_Rp = fun_Rp6, # from sA_pRp3
  fun_l = fun_l2, # from sA_lRl3
  fun_Rl = fun_Rl5, # from sA_lRl3
  fun_pl = detection3, # from sA_plRpl3
  fun_Rpl = fun_Rpl2, # from sA_plRpl3
  chunksize = 3e+7 
)
expect_equal(
  simulationsApply__8a[c(1,2,4,6,7)],
  simulationsApply__11[c(1,2,4,6,7)]
)
rm(simulationsApply__8a)
expect_equal(
  simulationsApply__11[[3]],
  fun_Rl5(simulationsApply__11[[2]])
)
expect_equal(
  simulationsApply__11[[5]],
  fun_Rp6(simulationsApply__11[[4]])
)
rm(simulationsApply__11)

## weight with default
fun_Rp7 = function(x, nout = 2, weight = 10){# how many of the models have Rsq < 0.5 (and slope <1)?
  c(mean(x[,3] <= 0.5, na.rm = TRUE) * weight[,2],
    mean(x[,2] <= 1, na.rm = TRUE))
}
fun_Rl6 = function(x, nout = 2, y = 2 * x, weight = 10){# extra parameter, with default: works
  c(0.5 * mean(y[,1], na.rm = TRUE),
    mean(x[,1], na.rm = TRUE)/ weight[,1])
}
simulationsApply__12 <- simulationsApply( # nout > 1, chunks
  simulations = radioactivePlumes_local,
  locations = locationsI,
  plumes = plumesI,
  fun = fun2, # from sA_fun3
  fun_p = fun_p2, # from sA_pRp3
  fun_Rp = fun_Rp7, # from sA_pRp3
  fun_l = fun_l2, # from sA_lRl3
  fun_Rl = fun_Rl6, # from sA_lRl3
  fun_pl = detection3, # from sA_plRpl3
  fun_Rpl = fun_Rpl2, # from sA_plRpl3
  chunksize = 3e+7 
)
expect_equal(
  simulationsApply__12[[3]],
  fun_Rl6(simulationsApply__12[[2]], weight = radioactivePlumes_local@locations@data[locationsI,, drop = FALSE])
)
expect_equal(
  simulationsApply__12[[5]],
  fun_Rp7(simulationsApply__12[[4]], weight = radioactivePlumes_local@plumes[plumesI,, drop = FALSE])
)
rm(simulationsApply__12)
## weight used but without default -> error
fun_Rp8 = function(x, nout = 2, weight){# how many of the models have Rsq < 0.5 (and slope <1)?
  c(mean(x[,3] <= 0.5, na.rm = TRUE) * weight[,2],
    mean(x[,2] <= 1, na.rm = TRUE))
}
fun_Rl7 = function(x, nout = 2, y = 2 * x, weight){# extra parameter, with default: works
  c(0.5 * mean(y[,1], na.rm = TRUE),
    mean(x[,1], na.rm = TRUE)/ weight[,2])
}
expect_error(
  simulationsApply( # nout > 1, chunks
    simulations = radioactivePlumes_local,
    locations = locationsI,
    plumes = plumesI,
    fun = fun2, # from sA_fun3
    fun_p = fun_p2, # from sA_pRp3
    fun_Rp = fun_Rp8, # from sA_pRp3
    fun_l = fun_l2, # from sA_lRl3
    fun_Rl = fun_Rl7, # from sA_lRl3
    fun_pl = detection3, # from sA_plRpl3
    fun_Rpl = fun_Rpl2, # from sA_plRpl3
    chunksize = 3e+7 
  )  
)

# test: generation an deleting of files
### if everything fits into memory, no files should be generated
### if result includes rasters not in memory, only these should be generated and kept
#### where are files saved? (working directory: s4p)
# generate a file
simulationsApply_pl1d_ <- simulationsApply( 
  simulations = radioactivePlumes_local,
  fun_pl = detection2,
  chunksize = 5e+6,
  nameSave = "tmp_simApply1",
  overwrite = TRUE
)
expect_true(
  file.exists(simulationsApply_pl1d_[["result_locationsplumes"]]@file@name)
)
expect_error(# file exists, not overwritten
  simulations = simulationsApply_pl1d_ <- simulationsApply( 
    radioactivePlumes_local,
    fun_pl = detection2,
    chunksize = 5e+6,
    nameSave = "tmp_simApply1",
    overwrite = FALSE
  )  
)
rm(simulationsApply_pl1d_); gc()

# generate and delete a file (nout of fun_pl wrong)
expect_warning(
  simulationsApply_plRpl1a_ <- simulationsApply( 
    radioactivePlumes_local,
    fun_pl = detection2b,
    nameSave = "tmp_simApply2",
    chunksize = 5e+6,
    overwrite = TRUE
  )  
)
expect_false(
  file.exists("tmp_simApply2.gri")
)
rm(simulationsApply_plRpl1a_)

# generate and delete a file (intermediate result because !isLocUnique)
simulationsApply__9_ = simulationsApply(
  simulations = radioactivePlumes_local,
  locations = locationsI,
  plumes = plumesI,
  fun_pl = detection2, # from sA_plRpl1
  chunksize = 5e+6,
  nameSave = "tmp_simApply3",
  overwrite = TRUE
)  
expect_false(
  file.exists("intermediateResult_locationsplumes.gri")  
)
rm(simulationsApply__9_)
})
