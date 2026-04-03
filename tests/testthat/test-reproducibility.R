context("reproducibility")

same_seg_site_haplo <- function(pop1, pop2, simParam) {
  hap1 <- unname(pullSegSiteHaplo(pop1, simParam = simParam))
  hap2 <- unname(pullSegSiteHaplo(pop2, simParam = simParam))
  isTRUE(all.equal(hap1, hap2, check.attributes = FALSE))
}

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

  test <- sum(hap1 - hap2) == 0
  # print(test)
  expect_true(test)

  test <- !sum(hap1 - hap3) == 0
  # print(test)
  expect_true(test)
})

test_that("crossing routines are reproducible in serial", {
  set.seed(101)
  founder <- quickHaplo(nInd = 4, nChr = 2, segSites = 20)
  SP <- SimParam$new(founder)
  pop <- newPop(founder, simParam = SP)
  SP$nThreads <- 1L

  crossPlan <- cbind(1:4, c(2, 3, 4, 1))

  set.seed(202)
  cross1 <- makeCross(pop, crossPlan, simParam = SP)
  set.seed(202)
  cross2 <- makeCross(pop, crossPlan, simParam = SP)
  expect_true(same_seg_site_haplo(cross1, cross2, SP))

  set.seed(303)
  dh1 <- makeDH(pop, nDH = 2, simParam = SP)
  set.seed(303)
  dh2 <- makeDH(pop, nDH = 2, simParam = SP)
  expect_true(same_seg_site_haplo(dh1, dh2, SP))

  set.seed(404)
  polyFounder <- quickHaplo(nInd = 3, nChr = 2, segSites = 20)
  SP2 <- SimParam$new(polyFounder)
  diploid <- newPop(polyFounder, simParam = SP2)
  tetraploid <- doubleGenome(diploid, simParam = SP2)
  SP2$quadProb <- 1 # to test the quadrivalent code path
  SP2$nThreads <- 1L

  set.seed(405)
  reduced1 <- reduceGenome(tetraploid, nProgeny = 2, simParam = SP2)
  set.seed(405)
  reduced2 <- reduceGenome(tetraploid, nProgeny = 2, simParam = SP2)
  expect_true(same_seg_site_haplo(reduced1, reduced2, SP2))
})

test_that("crossing routines are reproducible across OpenMP thread counts", {
  skip_if(getNumThreads() < 2)

  set.seed(501)
  founder <- quickHaplo(nInd = 4, nChr = 2, segSites = 20)
  SP <- SimParam$new(founder)
  pop <- newPop(founder, simParam = SP)

  crossPlan <- cbind(1:4, c(2, 3, 4, 1))

  SP$nThreads <- 1L
  set.seed(502)
  cross1 <- makeCross(pop, crossPlan, simParam = SP)
  SP$nThreads <- 2L
  set.seed(502)
  cross2 <- makeCross(pop, crossPlan, simParam = SP)
  expect_true(same_seg_site_haplo(cross1, cross2, SP))

  SP$nThreads <- 1L
  set.seed(503)
  dh1 <- makeDH(pop, nDH = 2, simParam = SP)
  SP$nThreads <- 2L
  set.seed(503)
  dh2 <- makeDH(pop, nDH = 2, simParam = SP)
  expect_true(same_seg_site_haplo(dh1, dh2, SP))

  set.seed(504)
  polyFounder <- quickHaplo(nInd = 3, nChr = 2, segSites = 20)
  SP2 <- SimParam$new(polyFounder)
  diploid <- newPop(polyFounder, simParam = SP2)
  tetraploid <- doubleGenome(diploid, simParam = SP2)
  SP2$quadProb <- 1 # to test the quadrivalent code path

  SP2$nThreads <- 1L
  set.seed(505)
  reduced1 <- reduceGenome(tetraploid, nProgeny = 2, simParam = SP2)
  SP2$nThreads <- 2L
  set.seed(505)
  reduced2 <- reduceGenome(tetraploid, nProgeny = 2, simParam = SP2)
  expect_true(same_seg_site_haplo(reduced1, reduced2, SP2))
})

test_that("standalone MapPop helpers respect explicit thread counts", {
  skip_if(getNumThreads() < 2)

  set.seed(601)
  founder <- quickHaplo(nInd = 4, nChr = 2, segSites = 6)

  segGeno1 <- pullSegSiteGeno(founder, nThreads = 1L)
  segGeno2 <- pullSegSiteGeno(founder, nThreads = 2L)
  expect_true(isTRUE(all.equal(unname(segGeno1), unname(segGeno2))))

  segHaplo1 <- pullSegSiteHaplo(founder, nThreads = 1L)
  segHaplo2 <- pullSegSiteHaplo(founder, nThreads = 2L)
  expect_true(isTRUE(all.equal(unname(segHaplo1), unname(segHaplo2))))

  markerNames <- c("1_1", "2_2")
  markerGeno1 <- pullMarkerGeno(founder, markers = markerNames, nThreads = 1L)
  markerGeno2 <- pullMarkerGeno(founder, markers = markerNames, nThreads = 2L)
  expect_true(isTRUE(all.equal(unname(markerGeno1), unname(markerGeno2))))

  markerHaplo1 <- pullMarkerHaplo(founder, markers = markerNames, nThreads = 1L)
  markerHaplo2 <- pullMarkerHaplo(founder, markers = markerNames, nThreads = 2L)
  expect_true(isTRUE(all.equal(unname(markerHaplo1), unname(markerHaplo2))))

  editedHaplo <- pullMarkerHaplo(founder, markers = "1_1", nThreads = 1L)
  editedHaplo[1, 1] <- 1L - editedHaplo[1, 1]
  founder1 <- setMarkerHaplo(founder, haplo = editedHaplo, nThreads = 1L)
  founder2 <- setMarkerHaplo(founder, haplo = editedHaplo, nThreads = 2L)
  edited1 <- pullMarkerHaplo(founder1, markers = "1_1", nThreads = 1L)
  edited2 <- pullMarkerHaplo(founder2, markers = "1_1", nThreads = 1L)
  expect_true(isTRUE(all.equal(unname(edited1), unname(edited2))))

  newHaplo <- matrix(
    rep(0:1, length.out = founder@nInd * founder@ploidy),
    ncol = 1
  )
  founder1 <- addSegSite(
    founder,
    siteName = "x",
    chr = 1,
    mapPos = 0.25,
    haplo = newHaplo,
    nThreads = 1L
  )
  founder2 <- addSegSite(
    founder,
    siteName = "x",
    chr = 1,
    mapPos = 0.25,
    haplo = newHaplo,
    nThreads = 2L
  )
  added1 <- pullMarkerHaplo(founder1, markers = "x", nThreads = 1L)
  added2 <- pullMarkerHaplo(founder2, markers = "x", nThreads = 1L)
  expect_true(isTRUE(all.equal(unname(added1), unname(added2))))
})
