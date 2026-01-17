context("reproducibility")

test_that("quickHaplo() and runMacs() are reproducible", {
  # ---- quickHaplo() test ----

  set.seed(123)
  pop <- quickHaplo(nInd = 4, nChr = 2, segSites = 10)
  hap1 <- pullSegSiteHaplo(pop)

  set.seed(123)
  pop <- quickHaplo(nInd = 4, nChr = 2, segSites = 10)
  hap2 <- pullSegSiteHaplo(pop)

  set.seed(124)
  pop <- quickHaplo(nInd = 4, nChr = 2, segSites = 10)
  hap3 <- pullSegSiteHaplo(pop)

  test <- sum(hap1 - hap2) == 0
  # print(test)
  expect_true(test)

  test <- !sum(hap1 - hap3) == 0
  # print(test)
  expect_true(test)

  # ---- runMacs(..., nThreads = 1) test ----

  set.seed(123)
  pop <- runMacs(
    nInd = 4,
    nChr = 2,
    segSites = 10,
    species = "GENERIC",
    nThreads = 1
  )
  hap1 <- pullSegSiteHaplo(pop)

  set.seed(123)
  pop <- runMacs(
    nInd = 4,
    nChr = 2,
    segSites = 10,
    species = "GENERIC",
    nThreads = 1
  )
  hap2 <- pullSegSiteHaplo(pop)

  set.seed(124)
  pop <- runMacs(
    nInd = 4,
    nChr = 2,
    segSites = 10,
    species = "GENERIC",
    nThreads = 1
  )
  hap3 <- pullSegSiteHaplo(pop)

  test <- sum(hap1 - hap2) == 0
  # print(test)
  expect_true(test)

  test <- !sum(hap1 - hap3) == 0
  # print(test)
  expect_true(test)

  # ---- runMacs(..., nThreads = 2) test ----

  set.seed(123)
  pop <- runMacs(
    nInd = 4,
    nChr = 2,
    segSites = 10,
    species = "GENERIC",
    nThreads = 2
  )
  hap1 <- pullSegSiteHaplo(pop)

  set.seed(123)
  pop <- runMacs(
    nInd = 4,
    nChr = 2,
    segSites = 10,
    species = "GENERIC",
    nThreads = 2
  )
  hap2 <- pullSegSiteHaplo(pop)

  set.seed(124)
  pop <- runMacs(
    nInd = 4,
    nChr = 2,
    segSites = 10,
    species = "GENERIC",
    nThreads = 2
  )
  hap3 <- pullSegSiteHaplo(pop)

  set.seed(123)
  pop <- runMacs(
    nInd = 4,
    nChr = 2,
    segSites = 10,
    species = "GENERIC",
    nThreads = NULL # so getNumThreads() kicks in
  )
  hap4 <- pullSegSiteHaplo(pop)

  test <- sum(hap1 - hap2) == 0
  # print(test)
  expect_true(test)

  test <- !sum(hap1 - hap3) == 0
  # print(test)
  expect_true(test)

  test <- sum(hap1 - hap4) == 0
  # print(test)
  expect_true(test)
})
