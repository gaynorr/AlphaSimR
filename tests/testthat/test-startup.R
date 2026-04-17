test_that("startup message reports OpenMP support", {
  msg = testthat::with_mocked_bindings(
    AlphaSimR:::.AlphaSimRStartupMessage(),
    getNumThreads = function() 4L,
    isOpenMPAvailable = function() TRUE,
    .package = "AlphaSimR"
  )

  expect_match(msg, "^AlphaSimR .* with OpenMP support")
  expect_match(msg, "using 4 threads by default")
  expect_match(msg, 'vignette\\("parallelization", package="AlphaSimR"\\)')
})

test_that("startup message handles singular thread wording", {
  msg = testthat::with_mocked_bindings(
    AlphaSimR:::.AlphaSimRStartupMessage(),
    getNumThreads = function() 1L,
    isOpenMPAvailable = function() TRUE,
    .package = "AlphaSimR"
  )

  expect_match(msg, "using 1 thread by default")
})

test_that("startup message reports missing OpenMP support", {
  msg = testthat::with_mocked_bindings(
    AlphaSimR:::.AlphaSimRStartupMessage(),
    getNumThreads = function() 1L,
    isOpenMPAvailable = function() FALSE,
    .package = "AlphaSimR"
  )

  expect_match(msg, "without OpenMP support")
  expect_match(msg, "single-threaded mode")
  expect_match(msg, 'vignette\\("parallelization", package="AlphaSimR"\\)')
})

test_that("OpenMP detection is consistent with getNumThreads", {
  hasOpenMP = AlphaSimR:::isOpenMPAvailable()

  expect_type(hasOpenMP, "logical")
  expect_length(hasOpenMP, 1)
  expect_true(getNumThreads() >= 1L)

  if (!hasOpenMP) {
    expect_equal(getNumThreads(), 1L)
  }
})
