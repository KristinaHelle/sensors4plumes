#################################################################
#             test replaceDefault                               #
#################################################################

test_that("replaceDefault",{
  # not a function
  fun1 = "a"
  fun2 = 5
  expect_warning(
    fun_1 <- replaceDefault(fun1)
  )
  expect_equal(
    fun1, fun_1[[1]]
  )
  expect_false(
    fun_1[[2]]
  )
  expect_warning(
    fun_2 <- replaceDefault(fun2, type = "fun.simulationsApply")
  )
  expect_equal(
    fun2, fun_2[[1]]
  )
  expect_false(
    fun_2[[2]]
  )
  # parameters missing
  ## required parameters
  fun3 = function(x){sum(x)}
  expect_warning(
    fun_3 <- replaceDefault(fun3, type = "fun.simulationsApply")
  )
  expect_equal(
    fun3, fun_3[[1]]
  )
  expect_false(
    fun_3[[2]]
  )
  ## parameters for which defaults are given
  fun4 = function(simulations, locations){1}
  expect_warning(
   fun_4 <- replaceDefault(fun4, type = "costFun", newDefaults = list(nout = 1, simulations = 10))  
  )
  expect_equal(
    fun_4[[1]], function(simulations = 10, locations){1} # fitting parameters replaced
  )  
  expect_false(
    fun_4[[2]]
  )
  expect_warning(
    fun_4a <- replaceDefault(fun4, type = "costFun", newDefaults = list(1, simulations = 10))  
  )
  expect_equal(
    fun_4, fun_4a
  )
  # defaults missing
  ## required
  fun5 = function(x = 3, nout){x}
  expect_warning(
    fun_5 <- replaceDefault(fun5, type = "fun.simulationsApply", newDefaults = list(x = 10))  
  )  
  expect_equal(
    fun_5[[1]], function(x = 10, nout){x} # default replaced
  )
  expect_false(
    fun_5[[2]]
  )
  # correct replacement (names match, no further tests)
  fun6 = function(x = 17, weight = "a", extraWeight = matrix(1:12, nrow = 3)){}
  fun_6 = replaceDefault(fun6, newDefaults = list(weight = 3, extraWeight = "c"), type =  "summaryFun.summaryPlumes") 
  expect_equal(
    fun_6[[1]],
    function(x = 17, weight = 3, extraWeight = "c"){}
  )
  expect_true(
    fun_6[[2]]
  )
})




# below old
# funName1 = function(locations, simulations, a = 2){}
# funName2 = 5
# funName3 = function(simulations, a = 2){}
# newDefaults1 = list(simulations = 5, locations = 10)
# newDefaults2 = list(simulations = 5, a = 10)
# newDefaults3 = list(a = 10, b = 20)
# newDefaults4 = list(simulations = 5, a = 10, b = 20)
# expect_error(
#   replaceDefault(fun = funName2, type  = "costFun", newDefaults = newDefaults2)
# )  
# #expect_warning(
# #  replaceDefault(fun = funName1, type  = "costFun", newDefaults = newDefaults3)
# #) # no warning as finName1 is already perfect; only b = 20 is ignored
# expect_warning(
#   replaceDefault(fun = funName3, type  = "costFun", newDefaults = newDefaults1)
# )
# funNameNew1 = replaceDefault(fun = funName1, type = "costFun", newDefaults = newDefaults1)
# funNameNew2 = replaceDefault(fun = funName1, type  = "costFun", newDefaults = newDefaults2)
# funNameNew3 = replaceDefault(fun = funName1, type  = "costFun", newDefaults = newDefaults4)
# expect_equal(
#   formals(funNameNew1),
#   pairlist(locations = 10, simulations = 5, a = 2)
# )
# #expect_equal(
# #  formals(funNameNew2),
# #  as.pairlist(alist(locations =, simulations = 5, a = 3))
# #)
# expect_equal(
#   formals(funNameNew1),
#   pairlist(locations = 10, simulations = 5, a = 2)
# )
# 
# rbga.bin = replaceDefault(
#   fun = rbga.bin, type = "rbga.bin",
#   newDefaults = list(size = 100, popSize = 50, zeroToOneRatio = 0.1, iters = 15))
# 
# 
