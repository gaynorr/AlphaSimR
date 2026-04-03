context("phenotypes")

test_that("asCategorical_converts_correctly",{
  cont = matrix(data = 0, nrow = 7, ncol = 3)
  cont[, 1] = c(-3, -2, -1, 0, 1, 2, 3)
  cont[, 2] = c(-3, -2, -1, 0, 1, 2, 3)
  cont[, 3] = c(-3, -2, -1, 0, 1, 2, 3)

  expect_equal(asCategorical(x = cont[, 1]),
               matrix(c(1, 1, 1, 2, 2, 2, 2)))
  expect_equal(asCategorical(x = cont[, 1], threshold = c(-1, 0, 1)),
               matrix(c(NA, NA, 1, 2, 2, NA, NA)))
  expect_equal(asCategorical(x = cont[, 1], threshold = c(-Inf, -1, 0, 1, Inf)),
               matrix(c(1, 1, 2, 3, 4, 4, 4)))

  expect_warning(asCategorical(x = cont[, 1], p = 0.5))
  expect_equal(suppressWarnings(asCategorical(x = cont[, 1], p = 0.5)),
               asCategorical(x = cont[, 1], p = c(0.5, 0.5)))

  trtMean = apply(X = cont, MARGIN = 2, FUN = mean)
  trtVar = apply(X = cont, MARGIN = 2, FUN = var)
  expect_equal(asCategorical(x = cont[, 1], p = c(0.5, 0.5), var = trtVar[1]),
               matrix(c(1, 1, 1, 2, 2, 2, 2)))
  expect_equal(asCategorical(x = cont[, 1], p = c(2/7, 1/7, 1/7, 3/7), var = trtVar[1]),
               matrix(c(1, 1, 2, 3, 4, 4, 4)))

  expect_error(asCategorical(x = cont))
  cont2 = asCategorical(x = cont,
                        threshold = list(NULL,
                                         c(-Inf, 0, Inf),
                                         NULL))
  cont2Exp = cont
  cont2Exp[, 2] = c(1, 1, 1, 2, 2, 2, 2)
  expect_equal(cont2, cont2Exp)

  expect_error(asCategorical(x = cont, p = c(0.5, 0.5)))
  pList = list(NULL, c(0.5, 0.5), NULL)
  expect_error(asCategorical(x = cont, p = pList))
  expect_error(asCategorical(x = cont, p = pList, mean = trtMean))
  cont2 = asCategorical(x = cont, p = pList, mean = trtMean, var=trtVar)
  cont2Exp = cont
  cont2Exp[, 2] = c(1, 1, 1, 2, 2, 2, 2)
  expect_equal(cont2, cont2Exp)
})
