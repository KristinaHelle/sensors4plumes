###########################################################################
# test cleanIndices                                                       #
###########################################################################

test_that("cleanIndices", {
  locTotal1 = 1:100
  locAll1 = c(0, 1, 1, sample.int(100, 60), 101, NA)
  locFix1 =  c(0, 2, 2, sample.int(100, 10), 101, NA)
  locInitial1 =  c(0, 3, 3, sample.int(100, 10), 101, NA)
  locInitial2 = matrix(c(0, 4, 4, sample(setdiff(1:100, locFix1), 10), 101, NA, locFix1[2],
                         0, 5, 5, sample(setdiff(1:100, locFix1), 10), 101, NA, locFix1[3],
                         0, 6, 6, sample(setdiff(1:100, locFix1), 10), 101, NA, locFix1[4]),
                     nrow = 3, byrow = TRUE)
  locInitial3 = matrix(c(0, 4, 4, sample(setdiff(1:100, locFix1), 10), 101, NA, locFix1[2],
                         0, 5, 5, sample(setdiff(1:100, locFix1), 10), 101, NA, locFix1[3],
                         0, 6, 6, sample(setdiff(1:90, locFix1), 10), 91, 92, 93),
                       nrow = 3, byrow = TRUE)
  locInitial4 = matrix(c(locFix1, 0,
                         locFix1, 0,
                         locFix1, 0),
                       nrow = 3, byrow = TRUE)
  expect_warning(
    locCleaned1 <- cleanIndices(locTotal1, locAll1, locFix1, locInitial1)
  )
  expect_warning(
    locCleaned2 <- cleanIndices(locTotal1, locAll1, locFix1, locInitial2)
  )
  expect_warning(
    locCleaned3 <- cleanIndices(locTotal1, locAll1, locFix1, locInitial3)
  )
  expect_warning(
    locCleaned4 <- cleanIndices(locTotal1, locAll1, locFix1, locInitial4)
  )
  # locInitial - vector
  expect_equal(
    locCleaned1$locationsAll,
    unique(setdiff(intersect(union(locAll1, locInitial1), locTotal1), locFix1))
  )
  expect_equal(
    locCleaned1$locationsFix,
    unique(intersect(locFix1, locTotal1))
  )
  expect_equal(
    locCleaned1$locationsInitial,
    unique(intersect(setdiff(locInitial1, locFix1), locTotal1))
  )
  
  # locInitial - matrix
  expect_equal(
    locCleaned2$locationsAll,
    unique(setdiff(intersect(union(locAll1, locInitial2), locTotal1), locFix1))
  )
  expect_equal(
    locCleaned3$locationsAll,
    unique(setdiff(intersect(union(locAll1, locInitial3), locTotal1), locFix1))
  )
  expect_equal(
    locCleaned4$locationsAll,
    unique(setdiff(intersect(union(locAll1, locInitial4), locTotal1), locFix1))
  )
  
  expect_equal(
    locCleaned2$locationsFix,
    unique(intersect(locFix1, locTotal1))
  )
  expect_equal(
    locCleaned3$locationsFix,
    locCleaned4$locationsFix
  )

#   expect_equal(
#     locCleaned2$locationsInitial[[3]],
#     unique(intersect(setdiff(locInitial2[[3]], locFix1), locTotal1)) 
#   )  
  i = sample.int(3,1)
  expect_equal(# same number of valid entries in each row -> all rows preserved
    locCleaned2$locationsInitial[i,],
    unique(intersect(setdiff(locInitial2[i,], locFix1), locTotal1)) 
  )  

  expect_equal(# varying numbers of valied entries per row -> not all rows preserved: knowledge from input - row3 becomes row1
    locCleaned3$locationsInitial[1,],
    unique(intersect(setdiff(locInitial3[3,], locFix1), locTotal1)) 
  ) 
  expect_equal(# no valid entries
    locCleaned4$locationsInitial,
    integer(0)
  )
})
