context("misc")

test_that("misc_and_miscPop",{
  founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)
  popOrig = newPop(founderPop, simParam=SP)
  multiPop = newMultiPop(popOrig, popOrig)

  expect_equal(popOrig@misc, list())
  expect_equal(popOrig@miscPop, list())

  popSub = popOrig[1]
  expect_equal(popSub@misc, list())
  expect_equal(popSub@miscPop, list())

  pop = popOrig
  pop@misc$vec = rnorm(n=2)
  pop@misc$mat = matrix(1:4, nrow=2)
  pop@misc$mtP = popOrig # hmm, should this actually be multiple pop objects or one with multiple individuals?
  pop@misc$mtLP = list(popOrig, popOrig)
  pop@misc$mtMP = multiPop
  # setting these miscPop elements just as an example - we are not testing them below,
  # because they get dropped in most/all operations on a pop
  pop@miscPop$vec = sum(pop@misc$vec)
  pop@miscPop$af = colMeans(pullSegSiteGeno(pop, simParam=SP))
  pop@miscPop$mat = matrix(1:4, nrow=2)

  popSub = pop[1]
  expect_equal(popSub@misc$vec, pop@misc$vec[1])
  expect_equal(popSub@misc$mat, pop@misc$mat[1, , drop=FALSE])
  expect_equal(popSub@misc$mtP, pop@misc$mtP[1])
  expect_equal(popSub@misc$mtLP, pop@misc$mtLP[1])
  expect_equal(popSub@misc$mtMP, pop@misc$mtMP[1])
  expect_equal(popSub@miscPop, list())

  popSub = pop[0]
  expect_equal(popSub@misc$vec, numeric(0))
  expect_equal(popSub@misc$mat, pop@misc$mat[0, , drop=FALSE])
  expect_equal(popSub@misc$mtP, newEmptyPop(simParam=SP))
  expect_equal(popSub@misc$mtLP, list())
  expect_equal(popSub@misc$mtMP, new("MultiPop", pops=list()))
  expect_equal(popSub@miscPop, list())

  popC = c(pop, pop)
  expect_equal(popC@misc$vec, c(pop@misc$vec, pop@misc$vec))
  expect_equal(popC@misc$mat, rbind(pop@misc$mat, pop@misc$mat))
  expect_equal(popC@misc$mtP, c(pop@misc$mtP, pop@misc$mtP))
  expect_equal(popC@misc$mtLP, c(pop@misc$mtLP, pop@misc$mtLP))
  expect_equal(popC@misc$mtMP, c(pop@misc$mtMP, pop@misc$mtMP))
  expect_equal(popC@miscPop, list())

  popC = c(pop, pop[1])
  expect_equal(popC@misc$vec, c(pop@misc$vec, pop@misc$vec[1]))
  expect_equal(popC@misc$mat, rbind(pop@misc$mat, pop@misc$mat[1, ]))
  expect_equal(popC@misc$mtP, c(pop@misc$mtP, pop@misc$mtP[1]))
  expect_equal(popC@misc$mtLP, c(pop@misc$mtLP, pop@misc$mtLP[1]))
  expect_equal(popC@misc$mtMP, c(pop@misc$mtMP, pop@misc$mtMP[1]))
  expect_equal(popC@miscPop, list())

  popA = popB = pop
  popA@misc = pop@misc[1]
  popB@misc = pop@misc[2]
  expect_warning(c(popA, popB)@misc,
                 regexp = "misc element names do not match - setting misc to an empty list!")

  popA@misc = pop@misc[1:2]
  popB@misc = pop@misc[2:1]
  expect_warning(c(popA, popB)@misc,
                 regexp = "misc element names do not match - setting misc to an empty list!")

  popA@misc = list()
  popB@misc = pop@misc[2:1]

  popList = list(popA, popB)
  expect_warning(c(popA, popB)@misc,
                 regexp = "number of misc elements differs - setting misc to an empty list!")

  popA@misc = pop@misc[2:1]
  popB@misc = list()
  expect_warning(c(popA, popB)@misc,
                 regexp = "number of misc elements differs - setting misc to an empty list!")
})
