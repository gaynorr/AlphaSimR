test_that("misc_and_miscPop",{
  founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)
  pop = newPop(founderPop, simParam=SP)

  popSub = pop[1]
  expect_equal(popSub@misc, list())
  expect_equal(popSub@miscPop, list())

  pop@misc$tmp1 = rnorm(n=2)
  pop@misc$tmp2 = rnorm(n=2)
  pop@miscPop$tmp1 = sum(pop@misc$tmp1)
  pop@miscPop$tmp2 = sum(pop@misc$tmp2)

  popSub = pop[1]
  expect_equal(popSub@misc$tmp1, pop@misc$tmp1[1])
  expect_equal(popSub@misc$tmp2, pop@misc$tmp2[1])
  expect_equal(popSub@miscPop, list())

  popSub = pop[0]
  expect_equal(popSub@misc$tmp1, numeric(0))
  expect_equal(popSub@misc$tmp2, numeric(0))
  expect_equal(popSub@miscPop, list())

  popC = c(pop, pop)
  expect_equal(popC@misc$tmp1, c(pop@misc$tmp1, pop@misc$tmp1))
  expect_equal(popC@misc$tmp2, c(pop@misc$tmp2, pop@misc$tmp2))
  expect_equal(popC@miscPop, list())
})
