# These tests will fail with a small probability and/or
# they are slower tests, so they are skipped on CRAN.
# If any of these test fail, rerun the tests multiple times.
# Frequent failures indicate a problem.
context("statistics")

test_that("addError", {
  skip_on_cran()
  gv = matrix(0, nrow = 10000, ncol = 2)
  varE = c(1, 1)
  pheno = AlphaSimR:::addError(gv = gv, varE = varE, reps = 1)
  expect_equal(var(pheno), diag(varE), tol = 0.1)
  varE = diag(2)
  pheno = AlphaSimR:::addError(gv = gv, varE = varE, reps = 1)
  expect_equal(var(pheno), varE, tol = 0.1)
  pheno = AlphaSimR:::addError(gv = gv, varE = varE, reps = 4)
  expect_equal(var(pheno), varE / 4, tol = 0.1)
  varE = 0.5 * diag(2) + 0.5
  pheno = AlphaSimR:::addError(gv = gv, varE = varE, reps = 1)
  expect_equal(var(pheno), varE, tol = 0.1)
  pheno = AlphaSimR:::addError(gv = gv, varE = varE, reps = 4)
  expect_equal(var(pheno), varE / 4, tol = 0.1)
  varE = 1.5 * diag(2) - 0.5
  pheno = AlphaSimR:::addError(gv = gv, varE = varE, reps = 1)
  expect_equal(var(pheno), varE, tol = 0.1)
  pheno = AlphaSimR:::addError(gv = gv, varE = varE, reps = 4)
  expect_equal(var(pheno), varE / 4, tol = 0.1)
  varE = matrix(c(1, 0.5, -0.5, 1), ncol = 2)
  expect_error(AlphaSimR:::addError(gv = gv, varE = varE, reps = 1))
})

test_that("sampleInt has correct invariants and approximately uniform output", {
  skip_on_cran()

  draws = AlphaSimR:::rngDiagnosticsSampleInt(
    n = 3L,
    N = 10L,
    reps = 1000L,
    seed = 101L
  )
  # Sampled indices should be sorted, unique, and stay within 0:(N - 1).
  expect_true(all(apply(draws, 1, function(x) identical(unname(sort(x)), x))))
  expect_true(all(apply(draws, 1, function(x) length(unique(x)) == length(x))))
  expect_true(min(draws) >= 0L)
  expect_true(max(draws) < 10L)

  empty = AlphaSimR:::rngDiagnosticsSampleInt(
    n = 0L,
    N = 4L,
    reps = 3L,
    seed = 102L
  )
  # This exercises the early return for n == 0.
  expect_equal(empty, matrix(integer(), nrow = 3L, ncol = 0L))

  full = AlphaSimR:::rngDiagnosticsSampleInt(
    n = 4L,
    N = 4L,
    reps = 3L,
    seed = 103L
  )
  # Sampling the full set should always return every value exactly once.
  expect_equal(full, matrix(rep(0:3, 3), nrow = 3, byrow = TRUE))

  # This exercises the explicit n > N error branch.
  expect_error(
    AlphaSimR:::rngDiagnosticsSampleInt(n = 5L, N = 4L, reps = 1L, seed = 104L),
    "n must be <= N"
  )

  one = AlphaSimR:::rngDiagnosticsSampleInt(
    n = 1L,
    N = 4L,
    reps = 40000L,
    seed = 105L
  )[, 1]
  probs = prop.table(table(factor(one, levels = 0:3)))
  # For n = 1, each of the four possible values should appear about 25% of the time.
  expect_true(all(abs(as.numeric(probs) - 0.25) <= 0.02))

  two = AlphaSimR:::rngDiagnosticsSampleInt(
    n = 2L,
    N = 4L,
    reps = 60000L,
    seed = 106L
  )
  observed = apply(two, 1, paste, collapse = "-")
  levels = apply(t(combn(0:3, 2)), 1, paste, collapse = "-")
  probs = prop.table(table(factor(observed, levels = levels)))
  # For n = 2, each of the six possible combinations should be about equally likely.
  expect_true(all(abs(as.numeric(probs) - 1 / 6) <= 0.02))
})

test_that("samplePoisson matches its target moments and zero probability", {
  skip_on_cran()

  lambdas = c(0.25, 1.0, 4.0)
  meanTol = c(0.02, 0.05, 0.15)
  zeroTol = c(0.01, 0.015, 0.02)

  for (i in seq_along(lambdas)) {
    lambda = lambdas[i]
    draws = AlphaSimR:::rngDiagnosticsSamplePoisson(
      lambda = lambda,
      reps = 100000L,
      seed = 101L + i
    )
    drawsR = rpois(
      n = 100000L,
      lambda = lambda
    )
    expect_equal(mean(draws), lambda, tolerance = meanTol[i])
    expect_equal(mean(drawsR), lambda, tolerance = meanTol[i])
    expect_equal(var(draws)[1, 1], lambda, tolerance = meanTol[i])
    expect_equal(var(drawsR), lambda, tolerance = meanTol[i])
    p <- ppois(q = 0L, lambda = lambda)
    expect_equal(mean(draws == 0L), p, tolerance = zeroTol[i])
    expect_equal(mean(drawsR == 0L), p, tolerance = zeroTol[i])
  }
})
