context("runMacs")

# MaCS draws the requested number of segregating sites while it simulates a
# chromosome, using reservoir sampling, rather than generating every site and
# subsampling afterwards. These tests cover the properties that behaviour has
# to keep. The statistical comparison of the retained sites against a uniform
# subset of all sites is too slow for the test suite and lives in
# macs_check.R at the top of the package.

test_that("runMacs returns exactly the requested number of sites", {
  set.seed(6001)
  pop = runMacs(nInd = 6, nChr = 3, segSites = 12, nThreads = 1)
  expect_equal(unname(sapply(pop@genMap, length)), rep(12L, 3L))
  expect_equal(unname(pop@nLoci), rep(12L, 3L))
})

test_that("the retained sites form a valid genetic map", {
  set.seed(6002)
  pop = runMacs(nInd = 6, nChr = 2, segSites = 15, nThreads = 1)
  for (chr in seq_along(pop@genMap)) {
    pos = pop@genMap[[chr]]
    # Sites come back in position order, starting at zero, with no repeats
    expect_false(is.unsorted(pos, strictly = TRUE))
    expect_equal(unname(pos[1]), 0)
    expect_false(any(duplicated(pos)))
  }
})

test_that("every retained site is segregating", {
  set.seed(6003)
  pop = runMacs(nInd = 8, nChr = 2, segSites = 15, nThreads = 1)
  SP = SimParam$new(pop)
  SP$nThreads = 1L
  basePop = newPop(pop, simParam = SP)
  freq = colMeans(pullSegSiteHaplo(basePop, simParam = SP))
  expect_true(all(freq > 0))
  expect_true(all(freq < 1))
})

test_that("asking for all sites returns more than a capped request", {
  set.seed(6004)
  capped = runMacs(nInd = 6, nChr = 1, segSites = 10, nThreads = 1)
  set.seed(6004)
  uncapped = runMacs(nInd = 6, nChr = 1, segSites = NULL, nThreads = 1)
  expect_equal(length(capped@genMap[[1]]), 10L)
  expect_gt(length(uncapped@genMap[[1]]), 10L)
})

test_that("requesting more sites than are generated is an error", {
  set.seed(6005)
  expect_error(runMacs(nInd = 4, nChr = 1, segSites = 1e6, nThreads = 1),
               "segSites")
})

test_that("runMacs site sampling is reproducible", {
  set.seed(6006)
  pop1 = runMacs(nInd = 6, nChr = 2, segSites = 12, nThreads = 1)
  set.seed(6006)
  pop2 = runMacs(nInd = 6, nChr = 2, segSites = 12, nThreads = 1)
  expect_equal(pop1@genMap, pop2@genMap)

  SP1 = SimParam$new(pop1); SP1$nThreads = 1L
  SP2 = SimParam$new(pop2); SP2$nThreads = 1L
  expect_equal(pullSegSiteHaplo(newPop(pop1, simParam = SP1), simParam = SP1),
               pullSegSiteHaplo(newPop(pop2, simParam = SP2), simParam = SP2))

  # A different seed gives a different map
  set.seed(6007)
  pop3 = runMacs(nInd = 6, nChr = 2, segSites = 12, nThreads = 1)
  expect_false(isTRUE(all.equal(pop1@genMap, pop3@genMap)))
})

test_that("runMacs site sampling does not depend on the thread count", {
  skip_on_cran()
  nThreads = getNumThreads()
  skip_if_not(nThreads > 1L, "only one thread available")

  set.seed(6008)
  pop1 = runMacs(nInd = 6, nChr = 4, segSites = 12, nThreads = 1)
  set.seed(6008)
  popN = runMacs(nInd = 6, nChr = 4, segSites = 12, nThreads = nThreads)
  expect_equal(pop1@genMap, popN@genMap)

  SP1 = SimParam$new(pop1); SP1$nThreads = 1L
  SPN = SimParam$new(popN); SPN$nThreads = 1L
  expect_equal(pullSegSiteHaplo(newPop(pop1, simParam = SP1), simParam = SP1),
               pullSegSiteHaplo(newPop(popN, simParam = SPN), simParam = SPN))
})
