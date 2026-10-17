test_that("mutate resolves unordered sites across chromosomes", {
  set.seed(42)
  founderPop = quickHaplo(nInd=2, nChr=3, segSites=4)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  pop = newPop(founderPop, simParam=SP)
  haplo = pullSegSiteHaplo(pop, simParam=SP)

  result = mutate(pop, mutRate=1, returnPos=TRUE, simParam=SP)
  expect_equal(pullSegSiteHaplo(result[[1]], simParam=SP), 1-haplo)
  expect_equal(pullSegSiteHaplo(pop, simParam=SP), haplo)
  positions = result[[2]]
  expect_equal(nrow(positions), 48L)
  expect_true(all(positions$chromosome %in% 1:3))
  expect_true(all(positions$site %in% 1:4))
  expect_equal(nrow(unique(positions)), 48L)

  unchanged = mutate(pop, mutRate=0, simParam=SP)
  expect_equal(pullSegSiteHaplo(unchanged, simParam=SP), haplo)
})
