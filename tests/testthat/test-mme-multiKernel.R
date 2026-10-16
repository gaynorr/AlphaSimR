context("mme-multiKernel")

# solveMKM and solveRRBLUPMK fit a single trait with several random terms by
# REML. There is no outside implementation to compare them against here, so
# these tests use things that must hold whatever the data: the solver's
# likelihood cannot be beaten by a brute force search, the estimates have to
# satisfy Henderson's mixed model equations at the variances returned, and
# the same model reached by different routes has to give the same answer.

# ---------------------------------------------------------------------------
# 1. Helpers
# ---------------------------------------------------------------------------

# Dense restricted log-likelihood, using the same convention as the solvers:
# -0.5*[(n-p)*log(2*pi) + log|V| + log|X'V^-1X| + y'Py], with p the rank of X
remlLL = function(vu, ve, y, X, Glist){
  y = c(y)
  n = length(y)
  V = diag(c(ve), n)
  for(i in seq_along(Glist)){
    V = V + c(vu)[i]*Glist[[i]]
  }
  R = chol(V)
  Vi = chol2inv(R)
  ViX = Vi %*% X
  XViX = crossprod(X, ViX)
  Rx = chol(XViX)
  P = Vi - ViX %*% chol2inv(Rx) %*% t(ViX)
  p = qr(X)$rank
  -0.5*((n-p)*log(2*pi) + 2*sum(log(diag(R))) + 2*sum(log(diag(Rx))) +
          drop(crossprod(y, P %*% y)))
}

# Maximizes remlLL over the log variances from several starting points. On
# the log scale no component can reach zero, so the maximum found is never
# above the constrained maximum the solver looks for.
bruteForceREML = function(y, X, Glist){
  k = length(Glist)
  vy = var(c(y))
  f = function(lt){
    v = exp(lt)
    val = tryCatch(remlLL(v[1:k], v[k+1], y, X, Glist),
                   error=function(e) -1e10)
    -val
  }
  best = list(value=Inf)
  for(s in c(0.1, 0.5, 0.9)){
    start = log(c(vy*s/(k*sapply(Glist, function(G) mean(diag(G)))),
                  vy*(1-s)))
    fit = optim(start, f, method="BFGS",
                control=list(maxit=1000, reltol=1e-12))
    fit = optim(fit$par, f, method="Nelder-Mead",
                control=list(maxit=5000, reltol=1e-14))
    if(fit$value < best$value){
      best = fit
    }
  }
  v = exp(best$par)
  list(LL=-best$value, Vu=v[1:k], Ve=v[k+1])
}

# Small simulated data set with two kernels and a signal in both
simTwoKernels = function(n=60, m1=40, m2=25, seed=9101){
  set.seed(seed)
  M1 = matrix(sample(-1:1, n*m1, replace=TRUE), n, m1)
  M2 = matrix(sample(0:1, n*m2, replace=TRUE), n, m2)
  y = 3 + M1 %*% rnorm(m1, sd=0.3) + M2 %*% rnorm(m2, sd=0.3) +
    rnorm(n, sd=1.5)
  list(y=y, X=matrix(1, n, 1), M1=M1, M2=M2, n=n)
}

# ---------------------------------------------------------------------------
# 2. Likelihood
# ---------------------------------------------------------------------------

test_that("OPTIMUM solveMKM is at least as good as a brute force REML search", {
  d = simTwoKernels()
  K1 = tcrossprod(d$M1)/ncol(d$M1) + diag(0.01, d$n)
  K2 = tcrossprod(d$M2)/ncol(d$M2) + diag(0.01, d$n)
  Z = diag(d$n)
  ans = solveMKM(d$y, d$X, list(Z, Z), list(K1, K2))
  ref = bruteForceREML(d$y, d$X, list(K1, K2))

  expect_gte(ans$LL, ref$LL - 1e-4)
  # LL is the convention above evaluated at the returned variances
  expect_equal(ans$LL, remlLL(ans$Vu, ans$Ve, d$y, d$X, list(K1, K2)),
               tolerance=1e-8)
})

test_that("OPTIMUM solveRRBLUPMK is at least as good as a brute force REML search", {
  d = simTwoKernels()
  G = list(tcrossprod(d$M1), tcrossprod(d$M2))
  ans = solveRRBLUPMK(d$y, d$X, list(d$M1, d$M2))
  ref = bruteForceREML(d$y, d$X, G)

  expect_gte(ans$LL, ref$LL - 1e-4)
  expect_equal(ans$LL, remlLL(ans$Vu, ans$Ve, d$y, d$X, G), tolerance=1e-8)
})

test_that("SINGLE a single kernel agrees with a single kernel REML search", {
  d = simTwoKernels(seed=9102)
  K = tcrossprod(d$M1)/ncol(d$M1) + diag(0.01, d$n)
  ans = solveMKM(d$y, d$X, list(diag(d$n)), list(K))
  ref = bruteForceREML(d$y, d$X, list(K))

  expect_gte(ans$LL, ref$LL - 1e-4)
  # The optimum is inside the parameter space here, so the two searches
  # have to find the same point
  expect_gt(c(ans$Vu), 0)
  expect_equal(c(ans$Vu), ref$Vu, tolerance=1e-3)
  expect_equal(c(ans$Ve), ref$Ve, tolerance=1e-3)

  rr = solveRRBLUPMK(d$y, d$X, list(d$M1))
  refRR = bruteForceREML(d$y, d$X, list(tcrossprod(d$M1)))
  expect_gte(rr$LL, refRR$LL - 1e-4)
  expect_equal(c(rr$Vu), refRR$Vu, tolerance=1e-3)
  expect_equal(c(rr$Ve), refRR$Ve, tolerance=1e-3)
})

# ---------------------------------------------------------------------------
# 3. Estimates
# ---------------------------------------------------------------------------

test_that("HENDERSON solveMKM estimates satisfy the mixed model equations", {
  # Records map onto 20 lines and 6 environments, so Z is not the identity
  set.seed(9103)
  n = 80
  line = sample(rep(1:20, 4))
  env = sample(rep(1:6, length.out=n))
  Z1 = outer(line, 1:20, "==")*1
  Z2 = outer(env, 1:6, "==")*1
  M = matrix(sample(-1:1, 20*30, replace=TRUE), 20, 30)
  K1 = tcrossprod(M)/30 + diag(0.05, 20)
  K2 = diag(6)
  X = cbind(1, rnorm(n))
  y = X %*% c(5, 1) + Z1 %*% (t(chol(K1)) %*% rnorm(20)) +
    Z2 %*% rnorm(6, sd=2) + rnorm(n)

  ans = solveMKM(y, X, list(Z1, Z2), list(K1, K2))
  expect_true(all(ans$Vu > 0))

  ve = c(ans$Ve)
  W = cbind(X, Z1, Z2)
  C = crossprod(W)
  i1 = ncol(X) + 1:20
  i2 = ncol(X) + 20 + 1:6
  C[i1, i1] = C[i1, i1] + ve/ans$Vu[1]*solve(K1)
  C[i2, i2] = C[i2, i2] + ve/ans$Vu[2]*solve(K2)
  sol = c(ans$beta, ans$u[[1]], ans$u[[2]])
  expect_equal(c(C %*% sol), c(crossprod(W, y)), tolerance=1e-8)
})

test_that("HENDERSON solveRRBLUPMK estimates satisfy the mixed model equations", {
  d = simTwoKernels(seed=9104)
  ans = solveRRBLUPMK(d$y, d$X, list(d$M1, d$M2))
  expect_true(all(ans$Vu > 0))

  ve = c(ans$Ve)
  W = cbind(d$X, d$M1, d$M2)
  C = crossprod(W)
  i1 = 1 + seq_len(ncol(d$M1))
  i2 = 1 + ncol(d$M1) + seq_len(ncol(d$M2))
  diag(C)[i1] = diag(C)[i1] + ve/ans$Vu[1]
  diag(C)[i2] = diag(C)[i2] + ve/ans$Vu[2]
  sol = c(ans$beta, ans$u[[1]], ans$u[[2]])
  expect_equal(c(C %*% sol), c(crossprod(W, d$y)), tolerance=1e-8)
})

test_that("ROUTES the marker form and the kernel form give the same model", {
  # M2 has more columns than there are records, which is the branch that
  # forms the kernel instead of working through the thin matrix
  set.seed(9105)
  n = 40
  M1 = matrix(sample(-1:1, n*15, replace=TRUE), n, 15)
  M2 = matrix(sample(0:1, n*70, replace=TRUE), n, 70)
  y = M1 %*% rnorm(15, sd=0.4) + M2 %*% rnorm(70, sd=0.15) + rnorm(n)
  X = matrix(1, n, 1)

  rr = solveRRBLUPMK(y, X, list(M1, M2))
  mk = solveMKM(y, X, list(diag(n), diag(n)),
                list(tcrossprod(M1), tcrossprod(M2)), tol=1e-8)

  expect_equal(c(rr$Vu), c(mk$Vu), tolerance=1e-4)
  expect_equal(c(rr$Ve), c(mk$Ve), tolerance=1e-4)
  expect_equal(rr$LL, mk$LL, tolerance=1e-6)
  # Genetic values from marker effects are the BLUPs of the kernel form
  expect_equal(c(M1 %*% rr$u[[1]]), c(mk$u[[1]]), tolerance=1e-3)
  expect_equal(c(M2 %*% rr$u[[2]]), c(mk$u[[2]]), tolerance=1e-3)
})

# ---------------------------------------------------------------------------
# 4. Shapes, boundaries and edge cases
# ---------------------------------------------------------------------------

test_that("SHAPES kernels with different column counts and a single column X", {
  set.seed(9106)
  n = 50
  Ms = list(matrix(rnorm(n*7), n), matrix(rnorm(n*30), n),
            matrix(rnorm(n*80), n))
  y = Ms[[1]] %*% rnorm(7, sd=0.5) + rnorm(n)
  X = matrix(1, n, 1)
  ans = solveRRBLUPMK(y, X, Ms)

  expect_equal(dim(ans$Vu), c(3L, 1L))
  expect_equal(dim(ans$Ve), c(1L, 1L))
  expect_equal(dim(ans$beta), c(1L, 1L))
  expect_equal(length(ans$u), 3L)
  expect_equal(sapply(ans$u, nrow), c(7L, 30L, 80L))
  expect_true(all(ans$Vu >= 0))
  expect_true(ans$Ve > 0)
  expect_true(all(is.finite(unlist(ans$u))))
})

test_that("BOUNDARY a component with no variance ends non-negative", {
  # The second term has no effect on y at all
  set.seed(9107)
  n = 60
  M1 = matrix(sample(-1:1, n*30, replace=TRUE), n)
  M2 = matrix(sample(0:1, n*30, replace=TRUE), n)
  X = matrix(1, n, 1)
  for(s in 1:5){
    y = M1 %*% rnorm(30, sd=0.3) + rnorm(n)
    rr = solveRRBLUPMK(y, X, list(M1, M2))
    mk = solveMKM(y, X, list(diag(n), diag(n)),
                  list(tcrossprod(M1)/30, tcrossprod(M2)/30))
    for(ans in list(rr, mk)){
      expect_true(all(ans$Vu >= 0))
      expect_true(all(is.finite(ans$Vu)))
      expect_true(is.finite(ans$LL))
      expect_false(anyNA(unlist(ans$u)))
      expect_false(anyNA(ans$beta))
    }
  }
})

test_that("BOUNDARY a kernel of zeros is held at zero", {
  d = simTwoKernels(seed=9108)
  ans = solveRRBLUPMK(d$y, d$X, list(d$M1, matrix(0, d$n, 5)))
  expect_equal(ans$Vu[2], 0)
  expect_equal(c(ans$u[[2]]), rep(0, 5))
  expect_gt(ans$Vu[1], 0)
})

test_that("RANK a dependent column of X changes nothing but its coefficient", {
  d = simTwoKernels(seed=9109)
  x2 = rnorm(d$n)
  X = cbind(1, x2)
  Xdup = cbind(1, x2, 2*x2 - 1)
  a = solveRRBLUPMK(d$y, X, list(d$M1, d$M2))
  b = solveRRBLUPMK(d$y, Xdup, list(d$M1, d$M2))
  expect_equal(b$LL, a$LL, tolerance=1e-8)
  expect_equal(c(b$Vu), c(a$Vu), tolerance=1e-6)
  expect_equal(c(Xdup %*% b$beta), c(X %*% a$beta), tolerance=1e-6)
})

test_that("RANK the intercept is kept wherever it is in X", {
  # The intercept is placed last, after two indicators that add up to it, so
  # the greedy choice would keep both indicators if it did not consider the
  # intercept first
  d = simTwoKernels(seed=9110)
  a = rep(c(1, 0), each=d$n/2)
  ref = solveRRBLUPMK(d$y, cbind(1, a), list(d$M1, d$M2))
  dep = solveRRBLUPMK(d$y, cbind(a, 1-a, 1), list(d$M1, d$M2))
  expect_equal(dep$LL, ref$LL, tolerance=1e-8)
  expect_equal(c(dep$Vu), c(ref$Vu), tolerance=1e-6)
  expect_equal(c(dep$beta), c(ref$beta[2], 0, ref$beta[1]), tolerance=1e-6)
})

test_that("MAXITER reaching the iteration limit gives a warning", {
  d = simTwoKernels(seed=9110)
  expect_warning(solveRRBLUPMK(d$y, d$X, list(d$M1, d$M2), maxIter=1L),
                 "did not converge")
  expect_warning(solveMKM(d$y, d$X, list(diag(d$n)),
                          list(tcrossprod(d$M1)), maxIter=1L),
                 "did not converge")
})

test_that("RANK a column of zeros in X is dropped", {
  d = simTwoKernels(seed=9112)
  a = solveRRBLUPMK(d$y, d$X, list(d$M1, d$M2))
  b = solveRRBLUPMK(d$y, cbind(d$X, 0), list(d$M1, d$M2))
  expect_equal(b$LL, a$LL, tolerance=1e-8)
  expect_equal(c(b$beta), c(a$beta, 0), tolerance=1e-6)
})

test_that("IDENTICAL two copies of one kernel fit the same model as one", {
  # The two components cannot be told apart, so only their sum is
  # determined and the average information matrix is singular. Whatever
  # split the solver settles on, the model it fits is the single kernel one.
  d = simTwoKernels(seed=9113)
  one = solveRRBLUPMK(d$y, d$X, list(d$M1))
  two = suppressWarnings(solveRRBLUPMK(d$y, d$X, list(d$M1, d$M1)))
  expect_true(all(two$Vu >= 0))
  expect_equal(two$LL, one$LL, tolerance=1e-6)
  expect_equal(sum(two$Vu), c(one$Vu), tolerance=1e-3)
  expect_equal(c(two$u[[1]] + two$u[[2]]), c(one$u[[1]]), tolerance=1e-3)
})

test_that("RESIDUAL a response that X fits exactly still gives an answer", {
  # A response of zeros leaves exactly no residual variance, so the
  # starting values need a fallback, and the components head for the
  # boundary
  d = simTwoKernels(seed=9114)
  y = matrix(0, d$n, 1)
  ans = suppressWarnings(solveRRBLUPMK(y, d$X, list(d$M1, d$M2)))
  expect_true(all(ans$Vu >= 0))
  expect_true(c(ans$Ve) > 0)
  expect_equal(c(ans$beta), 0, tolerance=1e-6)
  expect_false(anyNA(unlist(ans$u)))

  # A constant response leaves only rounding error behind, where the
  # likelihood is so flat that eventually no step improves it
  y = matrix(2, d$n, 1)
  ans = suppressWarnings(solveRRBLUPMK(y, d$X, list(d$M1, d$M2)))
  expect_true(all(ans$Vu >= 0))
  expect_true(c(ans$Ve) > 0)
  expect_equal(c(ans$beta), 2, tolerance=1e-6)
  expect_false(anyNA(unlist(ans$u)))
})

test_that("INPUT an indefinite kernel is rejected at the starting values", {
  # K has a positive diagonal but a large negative eigenvalue, so no
  # positive variance makes V positive definite
  d = simTwoKernels(seed=9115)
  K = diag(d$n) - (matrix(1, d$n, d$n) - diag(d$n))
  # The message depends on the side the solver chooses
  expect_error(solveMKM(d$y, d$X, list(diag(d$n)), list(K)),
               "not positive definite|must be positive semidefinite")
})

test_that("INPUT mismatched dimensions are rejected", {
  d = simTwoKernels(seed=9111)
  Z = diag(d$n)
  expect_error(solveRRBLUPMK(d$y, d$X, list(d$M1[-1, ])), "one row per record")
  expect_error(solveRRBLUPMK(cbind(d$y, d$y), d$X, list(d$M1)),
               "exactly one column")
  expect_error(solveRRBLUPMK(d$y, d$X[-1, , drop=FALSE], list(d$M1)),
               "same number of rows")
  expect_error(solveRRBLUPMK(d$y, d$X, list()), "at least one matrix")
  expect_error(solveRRBLUPMK(d$y, matrix(0, d$n, 1), list(d$M1)),
               "nonzero column")
  expect_error(solveRRBLUPMK(d$y, diag(d$n), list(d$M1)),
               "smaller than the number of records")
  expect_error(solveMKM(d$y, d$X, list(Z), list(diag(3))), "square")
  expect_error(solveMKM(d$y, d$X, list(Z), list()), "same length")
  expect_error(solveMKM(d$y, d$X, list(), list()), "at least one matrix")
  expect_error(solveMKM(d$y, d$X, list(Z[-1, ]), list(diag(d$n))),
               "one row per record")
})

test_that("INPUT missing values are rejected", {
  d = simTwoKernels(seed=9116)
  yNA = d$y
  yNA[3] = NA
  M1NA = d$M1
  M1NA[2, 2] = NA
  K = diag(d$n)
  K[1, 1] = NA
  expect_error(solveRRBLUPMK(yNA, d$X, list(d$M1)), "missing")
  expect_error(solveRRBLUPMK(d$y, d$X, list(M1NA)), "missing")
  expect_error(solveMKM(d$y, d$X, list(diag(d$n)), list(K)), "missing")
})

# ---------------------------------------------------------------------------
# 5. Record and coefficient sides
# ---------------------------------------------------------------------------

# The solvers work either with the variance matrix of the records or with
# Henderson's mixed model equations, whichever is cheaper. The two compute the
# same likelihood, scores and information, so on the same data they must take
# the same steps and reach the same estimates. The internal function
# .solveMultiKernel fits on a chosen side (1 records, 2 coefficients, 0
# automatic) and reports the side it used.
fitSide = function(y, X, Zlist, Klist, rrForm, side){
  AlphaSimR:::.solveMultiKernel(y, X, Zlist, Klist, rrForm, 200L, 1e-10,
                                as.integer(side))
}

expectSameFit = function(a, b){
  expect_equal(a$LL, b$LL, tolerance=1e-8)
  expect_equal(c(a$Vu), c(b$Vu), tolerance=1e-6)
  expect_equal(c(a$Ve), c(b$Ve), tolerance=1e-6)
  expect_equal(c(a$beta), c(b$beta), tolerance=1e-6)
  expect_equal(length(a$u), length(b$u))
  for(i in seq_along(a$u)){
    expect_equal(c(a$u[[i]]), c(b$u[[i]]), tolerance=1e-6)
  }
}

# Incidence matrix of n records spread over m levels, every level used
incidence = function(n, m){
  lev = c(seq_len(m), sample(m, n-m, replace=TRUE))
  Z = matrix(0, n, m)
  Z[cbind(seq_len(n), lev)] = 1
  Z
}

test_that("SIDES agree for the RR-BLUP form with fewer and more effects than records", {
  # 65 marker effects for 150 records, then for 40 records
  for(n in c(150, 40)){
    d = simTwoKernels(n=n, seed=9131)
    a = fitSide(d$y, d$X, list(d$M1, d$M2), list(), TRUE, 1)
    b = fitSide(d$y, d$X, list(d$M1, d$M2), list(), TRUE, 2)
    expect_equal(c(a$side, b$side), c(1, 2))
    expectSameFit(a, b)
  }
})

test_that("SIDES agree for a singular K and repeated records per level", {
  set.seed(9132)
  n = 90
  m = 20
  Z1 = incidence(n, m)
  Z2 = incidence(n, 15)
  B = matrix(rnorm(m*5), m, 5)
  K1 = tcrossprod(B)/5 # rank 5
  K2 = diag(15)
  y = Z1 %*% (B %*% rnorm(5)) + Z2 %*% rnorm(15, sd=0.5) + rnorm(n)
  X = cbind(1, rnorm(n))
  a = fitSide(y, X, list(Z1, Z2), list(K1, K2), FALSE, 1)
  b = fitSide(y, X, list(Z1, Z2), list(K1, K2), FALSE, 2)
  expectSameFit(a, b)
  # The BLUPs of a singular K lie in its column space, the span of B
  expect_equal(length(b$u[[1]]), m)
  expect_lt(max(abs(qr.resid(qr(B), c(b$u[[1]])))), 1e-8)
})

test_that("SIDES agree when a component is held at zero", {
  # The second term has no effect on y, so its variance often ends at zero,
  # and on the coefficient side it then leaves the equations
  set.seed(9133)
  n = 100
  M1 = matrix(sample(-1:1, n*30, replace=TRUE), n)
  M2 = matrix(sample(0:1, n*30, replace=TRUE), n)
  X = matrix(1, n, 1)
  nZero = 0
  for(s in 1:6){
    y = M1 %*% rnorm(30, sd=0.3) + rnorm(n)
    a = fitSide(y, X, list(M1, M2), list(), TRUE, 1)
    b = fitSide(y, X, list(M1, M2), list(), TRUE, 2)
    expectSameFit(a, b)
    if(b$Vu[2]==0){
      nZero = nZero + 1
      expect_equal(c(b$u[[2]]), rep(0, 30))
    }
  }
  expect_gt(nZero, 0)
})

test_that("SIDES agree with dependent columns in X", {
  d = simTwoKernels(n=150, seed=9134)
  a0 = rep(c(1, 0), each=75)
  X = cbind(a0, 1-a0, 1)
  a = fitSide(d$y, X, list(d$M1, d$M2), list(), TRUE, 1)
  b = fitSide(d$y, X, list(d$M1, d$M2), list(), TRUE, 2)
  expectSameFit(a, b)
  expect_equal(b$beta[2], 0)
})

test_that("OPTIMUM the coefficient side is at least as good as a brute force REML search", {
  d = simTwoKernels(n=200, seed=9135)
  ans = fitSide(d$y, d$X, list(d$M1, d$M2), list(), TRUE, 0)
  expect_equal(ans$side, 2)
  G = list(tcrossprod(d$M1), tcrossprod(d$M2))
  ref = bruteForceREML(d$y, d$X, G)
  expect_gte(ans$LL, ref$LL - 1e-4)
  expect_equal(ans$LL, remlLL(ans$Vu, ans$Ve, d$y, d$X, G), tolerance=1e-8)
})

test_that("AUTO the cheaper side is chosen", {
  # RR-BLUP form: 65 effects for 200 records, then for 40
  d = simTwoKernels(n=200, seed=9136)
  expect_equal(fitSide(d$y, d$X, list(d$M1, d$M2), list(), TRUE, 0)$side, 2)
  d = simTwoKernels(n=40, seed=9136)
  expect_equal(fitSide(d$y, d$X, list(d$M1, d$M2), list(), TRUE, 0)$side, 1)
  # Kernel form: 120 levels for 40 records, then 20 levels for 200 records
  set.seed(9137)
  n = 40
  m = 120
  M = matrix(rnorm(m*200), m)
  Z = cbind(diag(n), matrix(0, n, m-n))
  y = matrix(rnorm(n))
  expect_equal(fitSide(y, matrix(1, n, 1), list(Z),
                       list(tcrossprod(M)/200), FALSE, 0)$side, 1)
  n = 200
  y = matrix(rnorm(n))
  expect_equal(fitSide(y, matrix(1, n, 1), list(incidence(n, 20)),
                       list(diag(20)), FALSE, 0)$side, 2)
})

test_that("INPUT an indefinite K is rejected on the coefficient side", {
  set.seed(9138)
  n = 100
  m = 10
  K = diag(m) - (matrix(1, m, m) - diag(m))
  y = matrix(rnorm(n))
  expect_error(solveMKM(y, matrix(1, n, 1), list(incidence(n, m)), list(K)),
               "positive semidefinite")
  # The record side checks each K too, rather than only V
  expect_error(fitSide(y, matrix(1, n, 1), list(incidence(n, m)), list(K),
                       FALSE, 1), "positive semidefinite")
})

test_that("INPUT an indefinite K is rejected even when V stays positive definite", {
  # A small negative eigenvalue that the residual variance would hide in V
  set.seed(9140)
  n = 4
  y = matrix(rnorm(n))
  K = diag(c(1, 2, 3, -0.1))
  for(side in 1:2){
    expect_error(fitSide(y, matrix(1, n, 1), list(diag(n)), list(K), FALSE,
                         side), "positive semidefinite")
  }
})

test_that("SHIFT a large mean does not change the variance components", {
  # A fit on y+1e8 subtracts sums of squares of order 1e16 unless the mean is
  # taken off first. The variances must be unchanged and the intercept must
  # move by exactly the shift, on either side.
  d = simTwoKernels(n=120, seed=9141)
  for(side in 1:2){
    a = fitSide(d$y, d$X, list(d$M1, d$M2), list(), TRUE, side)
    b = fitSide(d$y+1e8, d$X, list(d$M1, d$M2), list(), TRUE, side)
    expect_equal(c(b$Vu), c(a$Vu), tolerance=1e-6)
    expect_equal(c(b$Ve), c(a$Ve), tolerance=1e-6)
    expect_equal(b$LL, a$LL, tolerance=1e-6)
    expect_equal(c(b$beta)-1e8, c(a$beta), tolerance=1e-6, scale=1)
  }
})

test_that("BOUNDARY no fit ends with a zero component that could raise the likelihood", {
  # Weak second terms often end at zero. Wherever they do, the likelihood
  # must be at least the best a search over positive values finds.
  skip_on_cran()
  set.seed(9142)
  n = 60
  M1 = matrix(sample(-1:1, n*20, replace=TRUE), n)
  M2 = matrix(sample(0:1, n*20, replace=TRUE), n)
  X = matrix(1, n, 1)
  G = list(tcrossprod(M1), tcrossprod(M2))
  y = M1 %*% rnorm(20, sd=0.3) + M2 %*% rnorm(20, sd=0.05) + rnorm(n)
  ref = bruteForceREML(y, X, G)
  for(side in 1:2){
    ans = fitSide(y, X, list(M1, M2), list(), TRUE, side)
    expect_gte(ans$LL, ref$LL - 1e-4)
  }
})

test_that("SIDES agree on a kernel of zeros, which has no effects to estimate", {
  set.seed(9139)
  n = 120
  Z = incidence(n, 15)
  y = Z %*% rnorm(15) + rnorm(n)
  Klist = list(diag(15), matrix(0, 15, 15))
  a = fitSide(y, matrix(1, n, 1), list(Z, Z), Klist, FALSE, 1)
  b = fitSide(y, matrix(1, n, 1), list(Z, Z), Klist, FALSE, 2)
  expectSameFit(a, b)
  expect_equal(b$Vu[2], 0)
  expect_equal(c(b$u[[2]]), rep(0, 15))
})

# ---------------------------------------------------------------------------
# 6. The genomic selection models that call solveRRBLUPMK
# ---------------------------------------------------------------------------

# A small diploid population with additive and dominance effects
mkPop = function(seed){
  set.seed(seed)
  founderPop = quickHaplo(nInd=60, nChr=2, segSites=20)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$restrSegSites(minQtlPerChr=10, minSnpPerChr=10, overlap=FALSE)
  SP$addTraitAD(nQtlPerChr=10, meanDD=0.5)
  SP$addSnpChip(nSnpPerChr=10)
  SP$setVarE(h2=0.5)
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)
  list(pop=pop, SP=SP)
}

test_that("RRBLUP_D matches solveRRBLUPMK fitted to the same matrices", {
  d = mkPop(9201)
  ans = RRBLUP_D(d$pop, simParam=d$SP)

  geno = pullSnpGeno(d$pop, simParam=d$SP)
  Ma = geno - 1
  Md = geno*(2 - geno)
  X = cbind(1, rowSums(Md))
  ref = solveRRBLUPMK(pheno(d$pop), X, list(Ma, Md))

  expect_equal(dim(ans@Vu), c(2L, 1L))
  expect_true(all(ans@Vu >= 0))
  expect_true(ans@Ve[1,1] > 0)
  expect_equal(c(ans@Vu), c(ref$Vu), tolerance=1e-8)
  expect_equal(c(ans@Ve), c(ref$Ve), tolerance=1e-8)
  expect_equal(c(ans@gv[[1]]@addEff), c(ref$u[[1]]), tolerance=1e-8)
  # Dominance effects of segregating markers carry the mean dominance
  # deviation, the coefficient on the inbreeding covariate
  seg = apply(geno, 2, function(x) length(unique(x)) > 1)
  expect_equal(c(ans@gv[[1]]@domEff)[seg],
               c(ref$u[[2]])[seg] + ref$beta[2], tolerance=1e-8)

  pred = setEBV(d$pop, ans, simParam=d$SP)
  expect_true(all(is.finite(c(ebv(pred)))))
})

test_that("RRBLUP_GCA matches solveRRBLUPMK fitted to the same matrices", {
  d = mkPop(9202)
  ans = RRBLUP_GCA(d$pop, simParam=d$SP)

  h1 = pullSnpHaplo(d$pop, haplo=1, simParam=d$SP)
  h2 = pullSnpHaplo(d$pop, haplo=2, simParam=d$SP)
  ref = solveRRBLUPMK(pheno(d$pop), matrix(1, nInd(d$pop), 1),
                      list(2*h1 - 1, 2*h2 - 1))

  expect_equal(dim(ans@Vu), c(2L, 1L))
  expect_true(all(ans@Vu >= 0))
  expect_equal(c(ans@Vu), c(ref$Vu), tolerance=1e-8)
  expect_equal(c(ans@Ve), c(ref$Ve), tolerance=1e-8)
  expect_equal(c(ans@female[[1]]@addEff), c(ref$u[[1]]), tolerance=1e-8)
  expect_equal(c(ans@male[[1]]@addEff), c(ref$u[[2]]), tolerance=1e-8)

  pred = setEBV(d$pop, ans, simParam=d$SP)
  expect_true(all(is.finite(c(ebv(pred)))))
})

test_that("RRBLUP_SCA matches solveRRBLUPMK fitted to the same matrices", {
  d = mkPop(9203)
  ans = RRBLUP_SCA(d$pop, simParam=d$SP)

  h1 = pullSnpHaplo(d$pop, haplo=1, simParam=d$SP)
  h2 = pullSnpHaplo(d$pop, haplo=2, simParam=d$SP)
  geno = h1 + h2
  Md = geno*(2 - geno)
  ref = solveRRBLUPMK(pheno(d$pop), cbind(1, rowSums(Md)),
                      list(2*h1 - 1, 2*h2 - 1, Md))

  expect_equal(dim(ans@Vu), c(3L, 1L))
  expect_true(all(ans@Vu >= 0))
  expect_true(ans@Ve[1,1] > 0)
  expect_equal(c(ans@Vu), c(ref$Vu), tolerance=1e-8)
  expect_equal(c(ans@Ve), c(ref$Ve), tolerance=1e-8)
  expect_equal(c(ans@gv[[1]]@addEff), c(ref$u[[1]]), tolerance=1e-8)
  expect_equal(c(ans@gv[[1]]@addEffMale), c(ref$u[[2]]), tolerance=1e-8)

  pred = setEBV(d$pop, ans, simParam=d$SP)
  expect_true(all(is.finite(c(ebv(pred)))))
})
