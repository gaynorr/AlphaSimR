context("phenotypes")

test_that("asLogNormal_converts_correctly", {
  x = matrix(data = c(-1, 0, 1, 0, 1, 2), nrow = 3, ncol = 2)

  expect_equal(asLogNormal(x = x[, 1]), matrix(exp(x[, 1])))
  expect_equal(
    asLogNormal(x = x[, 1], meanlog = 2),
    matrix(exp(2 + x[, 1]))
  )
  expect_equal(
    asLogNormal(x = x, meanlog = c(0, 3)),
    cbind(exp(0 + x[, 1]), exp(3 + x[, 2]))
  )
  expect_equal(
    asLogNormal(x = x, meanlog = list(NULL, 3)),
    cbind(x[, 1], exp(3 + x[, 2]))
  )

  expect_error(asLogNormal(x = x, meanlog = 0))
  expect_error(asLogNormal(x = x, meanlog = list(0)))
  expect_error(asLogNormal(x = x, meanlog = TRUE))
})

test_that("asCategorical_converts_correctly", {
  cont = matrix(data = 0, nrow = 7, ncol = 3)
  cont[, 1] = c(-3, -2, -1, 0, 1, 2, 3)
  cont[, 2] = c(-3, -2, -1, 0, 1, 2, 3)
  cont[, 3] = c(-3, -2, -1, 0, 1, 2, 3)

  expect_equal(asCategorical(x = cont[, 1]), matrix(c(1, 1, 1, 2, 2, 2, 2)))
  expect_equal(
    asCategorical(x = cont[, 1], threshold = c(-1, 0, 1)),
    matrix(c(NA, NA, 1, 2, 2, NA, NA))
  )
  expect_equal(
    asCategorical(x = cont[, 1], threshold = c(-Inf, -1, 0, 1, Inf)),
    matrix(c(1, 1, 2, 3, 4, 4, 4))
  )

  expect_warning(asCategorical(x = cont[, 1], p = 0.5))
  expect_equal(
    suppressWarnings(asCategorical(x = cont[, 1], p = 0.5)),
    asCategorical(x = cont[, 1], p = c(0.5, 0.5))
  )
  expect_error(asCategorical(x = cont[, 1], p = c(0.6, 0.6)))

  trtMean = apply(X = cont, MARGIN = 2, FUN = mean)
  trtVar = apply(X = cont, MARGIN = 2, FUN = var)
  expect_equal(
    asCategorical(x = cont[, 1], p = c(0.5, 0.5), var = trtVar[1]),
    matrix(c(1, 1, 1, 2, 2, 2, 2))
  )
  expect_equal(
    asCategorical(
      x = cont[, 1],
      p = c(2 / 7, 1 / 7, 1 / 7, 3 / 7),
      var = trtVar[1]
    ),
    matrix(c(1, 1, 2, 3, 4, 4, 4))
  )
  expect_equal(
    asCategorical(
      x = cont[, 1],
      threshold = c(-Inf, 0, Inf),
      include.lowest = TRUE,
      right = TRUE
    ),
    matrix(c(1, 1, 1, 1, 2, 2, 2))
  )

  expect_error(asCategorical(x = cont))
  cont2 = asCategorical(x = cont, threshold = list(NULL, c(-Inf, 0, Inf), NULL))
  cont2Exp = cont
  cont2Exp[, 2] = c(1, 1, 1, 2, 2, 2, 2)
  expect_equal(cont2, cont2Exp)

  expect_error(asCategorical(x = cont, p = c(0.5, 0.5)))
  pList = list(NULL, c(0.5, 0.5), NULL)
  expect_error(asCategorical(x = cont, p = pList))
  expect_error(asCategorical(x = cont, p = pList, mean = trtMean))
  cont2 = asCategorical(x = cont, p = pList, mean = trtMean, var = trtVar)
  cont2Exp = cont
  cont2Exp[, 2] = c(1, 1, 1, 2, 2, 2, 2)
  expect_equal(cont2, cont2Exp)
})

test_that("pop@gv and genParam(pop)@gv match", {
  # This test is here since we have two different code paths for these two
  # functionalities and we had one bug in one code path
  founderPop = quickHaplo(nInd = 10, nChr = 1, segSites = 10, ploidy = 1)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr = 10, mean = 0, var = 1, name = "addTraitA_allQTLs")
  SP$addTraitAE(
    nQtlPerChr = 10,
    relAA = 0,
    mean = 0,
    var = 1,
    useVarA = FALSE,
    name = "addTraitAE_allQTLs"
  )
  SP$addTraitA(nQtlPerChr = 2, mean = 0, var = 1, name = "addTraitA_2QTLs")
  SP$addTraitAE(
    nQtlPerChr = 2,
    relAA = 0,
    mean = 0,
    var = 1,
    useVarA = FALSE,
    name = "addTraitAE_2QTLs"
  )
  pop = newPop(founderPop, simParam = SP)
  diff = pop@gv - genParam(pop, simParam = SP)$gv
  test = abs(diff) < .Machine$double.eps^0.5
  expect_true(all(test))
})

makePhenoFinalizerTestSetup = function() {
  set.seed(101)
  founderPop = quickHaplo(nInd = 8, nChr = 1, segSites = 4)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr = 4, mean = c(0, 0), var = c(1, 2), corA = diag(2))
  SP$setVarE(varE = c(1, 1))

  list(founderPop = founderPop, SP = SP)
}

test_that("finalizePop can recode stored phenotypes", {
  setup = makePhenoFinalizerTestSetup()
  founderPop = setup$founderPop
  SP = setup$SP

  set.seed(201)
  pop = newPop(founderPop, simParam = SP)

  SP$finalizePop = function(pop, simParam = SP, targetTrait = 1L, ...) {
    pop@pheno[, targetTrait] = asCategorical(x = pheno(pop)[, targetTrait])
    pop
  }

  set.seed(201)
  finalizedPop = newPop(founderPop, simParam = SP)

  expected = pheno(pop)
  expected[, 1] = asCategorical(x = expected[, 1])

  expect_equal(pheno(finalizedPop), expected)

  set.seed(201)
  finalizedPop = newPop(founderPop, simParam = SP, targetTrait = 2L)

  expected = pheno(pop)
  expected[, 2] = asCategorical(x = expected[, 2])

  expect_equal(pheno(finalizedPop), expected)
})

test_that("finalizePheno can recode phenotypes during newPop", {
  setup = makePhenoFinalizerTestSetup()
  founderPop = setup$founderPop
  SP = setup$SP

  set.seed(202)
  pop = newPop(founderPop, simParam = SP)

  SP$finalizePheno = function(
    pheno,
    pop,
    simParam = SP,
    targetTrait = 1L,
    ...
  ) {
    pheno[, targetTrait] = asCategorical(x = pheno[, targetTrait])
    pheno
  }

  set.seed(202)
  finalizedPop = newPop(founderPop, simParam = SP)

  expected = pheno(pop)
  expected[, 1] = asCategorical(x = expected[, 1])

  expect_equal(pheno(finalizedPop), expected)

  set.seed(202)
  finalizedPop = newPop(founderPop, simParam = SP, targetTrait = 2L)

  expected = pheno(pop)
  expected[, 2] = asCategorical(x = expected[, 2])

  expect_equal(pheno(finalizedPop), expected)

  set.seed(202)
  finalizedPop = setPheno(finalizedPop, simParam = SP, targetTrait = 2L)

  expected = pheno(pop)
  expected[, 2] = asCategorical(x = expected[, 2])

  expect_equal(pheno(finalizedPop), expected)
})

test_that("finalizePheno can recode stored phenotypes via setPheno", {
  setup = makePhenoFinalizerTestSetup()
  founderPop = setup$founderPop
  SP = setup$SP
  pop = newPop(founderPop, simParam = SP)

  set.seed(203)
  expected = setPheno(pop, varE = c(1, 1), onlyPheno = TRUE, simParam = SP)

  SP$finalizePheno = function(
    pheno,
    pop,
    simParam = SP,
    targetTrait = 1L,
    ...
  ) {
    pheno[, targetTrait] = asCategorical(x = pheno[, targetTrait])
    pheno
  }

  set.seed(203)
  finalizedPop = setPheno(pop, varE = c(1, 1), simParam = SP, targetTrait = 2L)

  expected[, 2] = asCategorical(x = expected[, 2])

  expect_equal(pheno(finalizedPop), expected)
})

test_that("finalizePheno can recode onlyPheno output from setPheno", {
  setup = makePhenoFinalizerTestSetup()
  founderPop = setup$founderPop
  SP = setup$SP
  pop = newPop(founderPop, simParam = SP)

  set.seed(204)
  pheno = setPheno(pop, varE = c(1, 1), onlyPheno = TRUE, simParam = SP)

  SP$finalizePheno = function(
    pheno,
    pop,
    simParam = SP,
    targetTrait = 1L,
    ...
  ) {
    pheno[, targetTrait] = asCategorical(x = pheno[, targetTrait])
    pheno
  }

  set.seed(204)
  finalizedPheno = setPheno(
    pop,
    varE = c(1, 1),
    onlyPheno = TRUE,
    simParam = SP
  )

  expected = pheno
  expected[, 1] = asCategorical(x = expected[, 1])

  expect_equal(finalizedPheno, expected)

  set.seed(204)
  finalizedPheno = setPheno(
    pop,
    varE = c(1, 1),
    onlyPheno = TRUE,
    simParam = SP,
    targetTrait = 2L
  )

  expected = pheno
  expected[, 2] = asCategorical(x = expected[, 2])

  expect_equal(finalizedPheno, expected)
})

makeHybridPhenoFinalizerTestSetup = function() {
  setup = makePhenoFinalizerTestSetup()

  list(
    SP = setup$SP,
    pop = newPop(setup$founderPop[1:4], simParam = setup$SP),
    testers = newPop(setup$founderPop[5:8], simParam = setup$SP)
  )
}

calcExpectedSetPhenoGCA = function(
  pop,
  testers,
  simParam,
  targetTrait,
  varE = c(1, 1)
) {
  tmp = hybridCross(
    females = pop,
    males = testers,
    crossPlan = "testcross",
    simParam = simParam
  )
  y = setPheno(tmp, varE = varE, onlyPheno = TRUE, simParam = simParam)
  y[, targetTrait] = asCategorical(x = y[, targetTrait])

  female = factor(tmp@mother, levels = unique(tmp@mother))
  tmpAgg = aggregate(y ~ female, FUN = mean)
  unname(as.matrix(tmpAgg[, -1, drop = FALSE]))
}

calcExpectedSetPhenoProgTest = function(
  pop,
  testPop,
  simParam,
  targetTrait,
  nMatePerInd = 2L,
  varE = c(1, 1)
) {
  tmp = randCross2(
    females = pop,
    males = testPop,
    nCrosses = nInd(pop) * nMatePerInd,
    balance = TRUE,
    simParam = simParam
  )
  y = setPheno(tmp, varE = varE, onlyPheno = TRUE, simParam = simParam)
  y[, targetTrait] = asCategorical(x = y[, targetTrait])

  female = factor(tmp@mother, levels = pop@id)
  tmpAgg = aggregate(y ~ female, FUN = mean)
  unname(as.matrix(tmpAgg[, -1, drop = FALSE]))
}

test_that("finalizePheno can recode phenotypes via setPhenoGCA", {
  setup = makeHybridPhenoFinalizerTestSetup()
  SP = setup$SP
  pop = setup$pop
  testers = setup$testers

  set.seed(205)
  expected = calcExpectedSetPhenoGCA(
    pop = pop,
    testers = testers,
    simParam = SP,
    targetTrait = 2L
  )

  SP$finalizePheno = function(
    pheno,
    pop,
    simParam = SP,
    targetTrait = 1L,
    ...
  ) {
    pheno[, targetTrait] = asCategorical(x = pheno[, targetTrait])
    pheno
  }

  set.seed(205)
  finalizedPheno = setPhenoGCA(
    pop,
    testers,
    use = "pheno",
    varE = c(1, 1),
    onlyPheno = TRUE,
    simParam = SP,
    targetTrait = 2L
  )
  expect_equal(finalizedPheno, expected)

  set.seed(205)
  finalizedPop = setPhenoGCA(
    pop,
    testers,
    use = "pheno",
    varE = c(1, 1),
    simParam = SP,
    targetTrait = 2L
  )
  expect_equal(pheno(finalizedPop), expected)
})

test_that("finalizePheno can recode phenotypes via setPhenoProgTest", {
  setup = makeHybridPhenoFinalizerTestSetup()
  SP = setup$SP
  pop = setup$pop
  testers = setup$testers

  set.seed(206)
  expected = calcExpectedSetPhenoProgTest(
    pop = pop,
    testPop = testers,
    simParam = SP,
    targetTrait = 2L
  )

  SP$finalizePheno = function(
    pheno,
    pop,
    simParam = SP,
    targetTrait = 1L,
    ...
  ) {
    pheno[, targetTrait] = asCategorical(x = pheno[, targetTrait])
    pheno
  }

  set.seed(206)
  finalizedPheno = setPhenoProgTest(
    pop,
    testers,
    nMatePerInd = 2L,
    use = "pheno",
    varE = c(1, 1),
    onlyPheno = TRUE,
    simParam = SP,
    targetTrait = 2L
  )
  expect_equal(finalizedPheno, expected)

  set.seed(206)
  finalizedPop = setPhenoProgTest(
    pop,
    testers,
    nMatePerInd = 2L,
    use = "pheno",
    varE = c(1, 1),
    simParam = SP,
    targetTrait = 2L
  )
  expect_equal(pheno(finalizedPop), expected)
})
