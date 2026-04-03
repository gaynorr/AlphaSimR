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
