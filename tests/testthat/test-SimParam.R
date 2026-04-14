test_that("SimParam nThreads validates values and NULL resets to default", {
  founder <- quickHaplo(nInd = 2, nChr = 2, segSites = 4)
  SP <- SimParam$new(founder)
  SP$nThreads <- 1L
  pop <- newPop(founder, simParam = SP)
  expect_equal(SP$nThreads, 1L)

  SP$nThreads <- NULL
  expect_equal(SP$nThreads, getNumThreads())
  SP$nThreads <- 1L
  expect_silent(pullSegSiteGeno(pop, simParam = SP))

  expect_error(
    SP$nThreads <- 0,
    regexp = "single positive integer or NULL to reset"
  )
  expect_error(
    SP$nThreads <- 0L,
    regexp = "single positive integer or NULL to reset"
  )
  expect_error(
    SP$nThreads <- 1.5,
    regexp = "single positive integer or NULL to reset"
  )
  expect_error(
    SP$nThreads <- NA_integer_,
    regexp = "single positive integer or NULL to reset"
  )
})
