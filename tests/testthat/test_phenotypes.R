context("phenotypes")

test_that("asCategorical_converts_correctly",{
  cont = matrix(data = 0, nrow = 7, ncol = 3)
  cont[, 1] = c(-3, -2, -1, 0, 1, 2, 3)
  cont[, 2] = c(-3, -2, -1, 0, 1, 2, 3)
  cont[, 3] = c(-3, -2, -1, 0, 1, 2, 3)
  expect_equal(asCategorical(x = cont[, 2]),
               matrix(c(1, 1, 1, 2, 2, 2, 2)))
  expect_equal(asCategorical(x = cont[, 2], threshold = c(-1, 0, 1)),
               matrix(c(NA, NA, 1, 2, 2, NA, NA)))
  expect_equal(asCategorical(x = cont[, 2], threshold = c(-Inf, -1, 0, 1, Inf)),
               matrix(c(1, 1, 2, 3, 4, 4, 4)))
  expect_error(asCategorical(x = cont))
  cont2 = asCategorical(x = cont,
                        threshold = list(NULL,
                                         c(-Inf, 0, Inf),
                                         NULL))
  cont2Exp = cont
  cont2Exp[, 2] = c(1, 1, 1, 2, 2, 2, 2)
  expect_equal(cont2, cont2Exp)
})
