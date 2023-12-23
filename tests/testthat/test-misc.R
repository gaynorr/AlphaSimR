test_that("misc_and_miscPop",{
  founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)
  popOrig = newPop(founderPop, simParam=SP)
  multiPop = newMultiPop(popOrig, popOrig)

  popSub = popOrig[1]
  expect_equal(popSub@misc, list())
  expect_equal(popSub@miscPop, list())

  pop = popOrig
  pop@misc$vec = rnorm(n=2)
  pop@misc$mtP = popOrig # hmm, should this actually be multiple pop objects or one with multiple individuals?
  pop@misc$mtLP = list(popOrig, popOrig)
  pop@misc$mtMP = multiPop
  pop@miscPop$vec = sum(pop@misc$vec)
  pop@miscPop$af = colMeans(pullSegSiteGeno(pop))

  popSub = pop[1]
  expect_equal(popSub@misc$vec, pop@misc$vec[1])
  expect_equal(popSub@misc$mtP, pop@misc$mtP[1])
  expect_equal(popSub@misc$mtLP, pop@misc$mtLP[1])
  expect_equal(popSub@misc$mtMP, pop@misc$mtMP[1])
  expect_equal(popSub@miscPop, list())

  popSub = pop[0]
  expect_equal(popSub@misc$vec, numeric(0))
  expect_equal(popSub@misc$mtP, newEmptyPop())
  expect_equal(popSub@misc$mtLP, list())
  expect_equal(popSub@misc$mtMP, ???) # what should we test here against?
  expect_equal(popSub@miscPop, list())

  popC = c(pop, pop) # this now breaks with validity check
  # Error in validObject(.Object) :
  #   invalid class “Pop” object: any(nInd!=sapply(misc, length))
  # --> because length(pop@misc$vec) is 2, but length(pop@misc$mtP) is 1

  pop@misc$mtP = NULL # removing mtP for now to test how it works with mtLP and mtMP
  popC = c(pop, pop) # this now breaks with validity check
  # Error in validObject(.Object) :
  #   invalid class “Pop” object: any(nInd!=sapply(misc, length))
  # --> because length(pop@misc$vec) is 2, but length(pop@misc$mtMP) is 1
  # the later should really be 2!

  pop@misc$mtMP = NULL # removing mtMP for now to test how it works with mtLP
  popC = c(pop, pop) # this now breaks with validity check

  expect_equal(popC@misc$vec, c(pop@misc$vec, pop@misc$vec))
  expect_equal(popC@misc$mtP, c(pop@misc$mtP, pop@misc$mtP)) # failing atm
  expect_equal(popC@misc$mtLP, c(pop@misc$mtLP, pop@misc$mtLP)) # failing atm
  expect_equal(popC@misc$mtMP, c(pop@misc$mtMP, pop@misc$mtMP)) # failing atm
  expect_equal(popC@miscPop, list())
})
