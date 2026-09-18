context("miscMath")

# selInt, smithHazel, selIndex and usefulness are small pieces of arithmetic
# with answers that can be worked out independently. None of them had a test.
# They are the cheapest things in the package to check and among the easiest
# to get subtly wrong, because a transposed matrix or a misplaced inverse
# still returns a plausible looking number.

test_that("selInt matches the standard normal selection intensity", {
  # i = phi(z) / p where z is the truncation point leaving proportion p
  for(p in c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.9)){
    z = qnorm(1 - p)
    expect_equal(selInt(p), dnorm(z) / p, tolerance=1e-10,
                 info=paste("p =", p))
  }

  # Values from Falconer and Mackay's table of selection intensities
  expect_equal(selInt(0.01), 2.665, tolerance=1e-3)
  expect_equal(selInt(0.05), 2.063, tolerance=1e-3)
  expect_equal(selInt(0.10), 1.755, tolerance=1e-3)
  expect_equal(selInt(0.20), 1.400, tolerance=1e-3)
  expect_equal(selInt(0.50), 0.798, tolerance=1e-3)

  # Selecting everyone gives the mean of the distribution, which is zero
  expect_equal(selInt(1), 0, tolerance=1e-10)

  # Intensity rises as the selected proportion falls
  p = c(0.5, 0.2, 0.1, 0.05, 0.01)
  expect_true(all(diff(selInt(p)) > 0))

  # Vectorized, because it is built from vectorized pieces
  expect_equal(selInt(p), sapply(p, selInt), tolerance=1e-12)
})

test_that("smithHazel solves the index equations", {
  G = 1.5*diag(2) - 0.5
  E = diag(2)
  P = G + E
  wt = c(1, 1)

  b = smithHazel(wt, G, P)

  # The defining property: P b = G w
  expect_equal(c(P %*% b), c(G %*% wt), tolerance=1e-10)

  # Shape of the answer
  expect_true(is.matrix(b))
  expect_equal(nrow(b), 2L)
  expect_equal(ncol(b), 1L)

  # Two uncorrelated traits of equal heritability get equal weight when the
  # economic weights are equal, and weight in proportion when they are not
  G0 = diag(2)
  P0 = 2*diag(2)
  expect_equal(c(smithHazel(c(1,1), G0, P0)), c(0.5, 0.5), tolerance=1e-10)
  expect_equal(c(smithHazel(c(2,1), G0, P0)), c(1.0, 0.5), tolerance=1e-10)

  # A trait with no genetic variance gets no weight
  Gz = diag(c(1, 0))
  Pz = diag(c(2, 2))
  expect_equal(c(smithHazel(c(1,1), Gz, Pz))[2], 0, tolerance=1e-10)

  # Scaling the economic weights scales the index weights
  expect_equal(c(smithHazel(2*wt, G, P)), 2*c(smithHazel(wt, G, P)),
               tolerance=1e-10)

  # With one trait it reduces to the ratio of genetic to phenotypic variance,
  # which is the heritability
  expect_equal(c(smithHazel(1, matrix(0.5), matrix(2))), 0.25,
               tolerance=1e-10)
})

test_that("selIndex forms the weighted sum of trait values", {
  Y = matrix(c(1, 2, 3,
               4, 5, 6), nrow=3)
  b = c(2, 10)

  expect_equal(c(selIndex(Y, b)), c(1*2 + 4*10, 2*2 + 5*10, 3*2 + 6*10))
  expect_equal(nrow(selIndex(Y, b)), 3L)
  expect_equal(ncol(selIndex(Y, b)), 1L)

  # Scaling centres and standardizes each column first, so the result has
  # mean zero
  scaled = selIndex(Y, b, scale=TRUE)
  expect_equal(mean(scaled), 0, tolerance=1e-10)
  expect_equal(c(scaled), c(scale(Y) %*% b), tolerance=1e-10)

  # Weighting only the first trait ignores the second entirely
  expect_equal(c(selIndex(Y, c(1, 0))), Y[,1])
})

test_that("selIndex and smithHazel work together as a selection criterion", {
  set.seed(9301)
  founderPop = quickHaplo(nInd=40, nChr=2, segSites=20)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  G = 1.5*diag(2) - 0.5
  SP$addTraitA(10, mean=c(0,0), var=c(1,1), corA=G)
  SP$setVarE(h2=c(0.5,0.5))
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)

  b = smithHazel(c(1,1), varG(pop), varP(pop))
  expect_true(all(is.finite(c(b))))

  selected = selectInd(pop, nInd=10, trait=selIndex, simParam=SP, b=b)
  expect_equal(nInd(selected), 10L)
  expect_true(all(selected@id %in% pop@id))

  # The selected individuals are the ten with the highest index values
  index = c(selIndex(pheno(pop), b))
  best = pop@id[order(index, decreasing=TRUE)[1:10]]
  expect_true(setequal(selected@id, best))
})

test_that("usefulness is the mean of the selected tail", {
  set.seed(9401)
  founderPop = quickHaplo(nInd=50, nChr=1, segSites=20)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10)
  SP$setVarE(h2=0.5)
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)

  g = c(gv(pop))

  # Worked out by hand: sort, keep the top ceiling(p*n), average
  for(p in c(0.1, 0.25, 0.5)){
    n = ceiling(p * length(g))
    byHand = mean(sort(g, decreasing=TRUE)[1:n])
    expect_equal(usefulness(pop, use="gv", p=p, simParam=SP), byHand,
                 tolerance=1e-10, info=paste("p =", p))
  }

  # Selecting downward takes the other tail
  n = ceiling(0.1 * length(g))
  expect_equal(usefulness(pop, use="gv", p=0.1, selectTop=FALSE, simParam=SP),
               mean(sort(g, decreasing=FALSE)[1:n]), tolerance=1e-10)

  # Taking everyone gives the population mean
  expect_equal(usefulness(pop, use="gv", p=1, simParam=SP), mean(g),
               tolerance=1e-10)

  # A smaller selected proportion cannot do worse than a larger one
  expect_gte(usefulness(pop, use="gv", p=0.1, simParam=SP),
             usefulness(pop, use="gv", p=0.5, simParam=SP))

  # It reads genetic values by default
  expect_equal(usefulness(pop, p=0.5, simParam=SP),
               usefulness(pop, use="gv", p=0.5, simParam=SP))
  expect_false(isTRUE(all.equal(
    usefulness(pop, use="gv", p=0.5, simParam=SP),
    usefulness(pop, use="pheno", p=0.5, simParam=SP))))
})

test_that("attrition removes individuals at about the rate asked for", {
  set.seed(9501)
  founderPop = quickHaplo(nInd=500, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(5)
  pop = newPop(founderPop, simParam=SP)

  # Nothing lost and everything lost are exact
  expect_equal(nInd(attrition(pop, p=0)), nInd(pop))
  expect_equal(nInd(attrition(pop, p=1)), 0L)

  # Survivors are a subset of the original, in the original order
  kept = attrition(pop, p=0.5)
  expect_true(all(kept@id %in% pop@id))
  expect_false(any(duplicated(kept@id)))
  expect_equal(kept@id, pop@id[pop@id %in% kept@id])

  # The loss rate is binomial, so 500 individuals put the count well inside
  # this range and the test does not depend on the seed holding exactly
  set.seed(9502)
  n = nInd(attrition(pop, p=0.5))
  expect_gt(n, 200L)
  expect_lt(n, 300L)
})

test_that("getPed returns the pedigree of a population", {
  set.seed(9601)
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackPed(TRUE)
  SP$addTraitA(5)
  parents = newPop(founderPop, simParam=SP)
  progeny = randCross(parents, nCrosses=5, simParam=SP)

  ped = getPed(progeny)
  expect_equal(nrow(ped), nInd(progeny))
  expect_true(all(c("id","mother","father") %in% colnames(ped)))

  # Every parent named is one of the individuals that were crossed
  expect_true(all(as.character(ped$mother) %in% parents@id))
  expect_true(all(as.character(ped$father) %in% parents@id))
  expect_equal(as.character(ped$id), progeny@id)
})
