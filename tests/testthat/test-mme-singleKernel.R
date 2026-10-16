context("mme-singleKernel")

# solveUVM fits y = X*beta + Z*u + e with u ~ N(0, K*Vu) by REML. The tests
# check it against two things that must hold whatever the implementation: a
# brute force REML fit written here from the textbook likelihood, and
# Henderson's mixed model equations at the returned variances. The cases are
# chosen to reach both of its decompositions, the one on Z*K*Z' when there
# are at least as many levels as records and the one through K when there are
# fewer, together with a singular K and a single fixed effect.

# 1. Helpers --------------------------------------------------------------

# Dense restricted log-likelihood, excluding the log|X'X| term, at variances
# Vu and Ve. This is the formula solveUVM documents for LL.
denseREML = function(Vu, Ve, y, X, Z, K){
  n = nrow(X)
  p = ncol(X)
  V = Vu*Z%*%K%*%t(Z) + diag(Ve, n)
  Vi = solve(V)
  XVX = t(X)%*%Vi%*%X
  P = Vi - Vi%*%X%*%solve(XVX, t(X)%*%Vi)
  ll = (n-p)*log(2*pi) +
    c(determinant(V)$modulus) +
    c(determinant(XVX)$modulus) +
    c(t(y)%*%P%*%y)
  return(-0.5*ll)
}

# Maximizes denseREML over the log variances with a general purpose
# optimizer, starting from an even split of the phenotypic variance so that
# nothing is borrowed from solveUVM
bruteREML = function(y, X, Z, K){
  start = rep(log(var(c(y))/2), 2)
  obj = function(par) denseREML(exp(par[1]), exp(par[2]), y, X, Z, K)
  fit = optim(start, obj, method="BFGS",
              control=list(fnscale=-1, reltol=1e-14, maxit=1000))
  fit = optim(fit$par, obj, method="Nelder-Mead",
              control=list(fnscale=-1, reltol=1e-14, maxit=5000))
  return(list(Vu=exp(fit$par[1]), Ve=exp(fit$par[2]), LL=fit$value))
}

# Simulates a data set from the model
simUVM = function(n, X, Z, K, Vu=2, Ve=1, seed=1){
  set.seed(seed)
  m = ncol(Z)
  e = eigen(K, symmetric=TRUE)
  lambda = pmax(e$values, 0)
  u = e$vectors%*%(sqrt(lambda*Vu)*rnorm(m))
  beta = seq_len(ncol(X))
  y = X%*%beta + Z%*%u + rnorm(n, sd=sqrt(Ve))
  return(matrix(y, ncol=1))
}

# A random positive definite kernel, or one of reduced rank when rank < m
makeK = function(m, rank=m, seed=2){
  set.seed(seed)
  A = matrix(rnorm(m*rank), m, rank)
  K = tcrossprod(A)/rank
  if(rank==m){
    K = K + diag(0.1, m)
  }
  return(K)
}

# Incidence matrix giving each of n records one of m levels, every level used
makeZ = function(n, m, seed=3){
  set.seed(seed)
  lev = sample(rep_len(seq_len(m), n))
  Z = matrix(0, n, m)
  Z[cbind(seq_len(n), lev)] = 1
  return(Z)
}

# Checks the fit against the brute force reference and against Henderson's
# equations. The equations are used in the form multiplied through by K,
#
#   X'X*beta + X'Z*u              = X'y
#   K*Z'X*beta + (K*Z'Z + delta*I)*u = K*Z'y
#
# for delta = Ve/Vu, which holds whether or not K can be inverted. When K is
# singular the BLUPs must also lie in its column space.
checkUVM = function(y, X, Z, K){
  ans = solveUVM(y, X, Z, K)
  ref = bruteREML(y, X, Z, K)

  # The likelihood reported must be the likelihood at the reported variances
  expect_equal(ans$LL, denseREML(ans$Vu, ans$Ve, y, X, Z, K),
               tolerance=1e-8, scale=1)
  # and must be the maximum, up to how precisely each search can locate it
  expect_gt(ans$LL, ref$LL - 1e-6)
  expect_equal(ans$LL, ref$LL, tolerance=1e-6, scale=1)
  expect_equal(ans$Vu, ref$Vu, tolerance=1e-3)
  expect_equal(ans$Ve, ref$Ve, tolerance=1e-3)

  beta = ans$beta
  u = ans$u
  delta = ans$Ve/ans$Vu
  expect_equal(dim(beta), c(ncol(X), 1L))
  expect_equal(dim(u), c(ncol(Z), 1L))
  lhs1 = t(X)%*%X%*%beta + t(X)%*%Z%*%u
  lhs2 = K%*%t(Z)%*%X%*%beta + (K%*%t(Z)%*%Z + diag(delta, ncol(Z)))%*%u
  expect_equal(c(lhs1), c(t(X)%*%y), tolerance=1e-8)
  expect_equal(c(lhs2), c(K%*%t(Z)%*%y), tolerance=1e-8)

  e = eigen(K, symmetric=TRUE)
  nullK = e$vectors[, e$values < 1e-8*max(e$values), drop=FALSE]
  if(ncol(nullK) > 0){
    expect_equal(c(t(nullK)%*%u), rep(0, ncol(nullK)),
                 tolerance=1e-8, scale=max(abs(u)))
  }
  invisible(ans)
}

# 2. solveUVM -------------------------------------------------------------

test_that("FIT one record per level with a full rank K", {
  set.seed(13)
  n = 30
  X = cbind(1, rnorm(n))
  Z = diag(n)
  K = makeK(n)
  y = simUVM(n, X, Z, K)
  checkUVM(y, X, Z, K)
})

test_that("FIT one record per level with a singular K", {
  n = 30
  X = cbind(1, seq_len(n)/n)
  Z = diag(n)
  K = makeK(n, rank=8)
  y = simUVM(n, X, Z, K, seed=4)
  checkUVM(y, X, Z, K)
})

test_that("FIT repeated records per level with a single fixed effect", {
  # 12 levels and 48 records, so the decomposition goes through K
  n = 48
  m = 12
  X = matrix(1, n, 1)
  Z = makeZ(n, m)
  K = makeK(m)
  y = simUVM(n, X, Z, K, seed=5)
  checkUVM(y, X, Z, K)
})

test_that("FIT repeated records per level with a singular K", {
  # K has rank 6 of 15, so some of its eigenvectors are discarded before
  # Z*K*Z' is decomposed
  set.seed(14)
  n = 45
  m = 15
  X = cbind(1, rnorm(n))
  Z = makeZ(n, m, seed=6)
  K = makeK(m, rank=6)
  y = simUVM(n, X, Z, K, seed=7)
  checkUVM(y, X, Z, K)
})

test_that("FIT more levels than records", {
  # Levels with no records still receive BLUPs through K
  n = 20
  m = 30
  X = matrix(1, n, 1)
  Z = matrix(0, n, m)
  Z[cbind(seq_len(n), seq_len(n))] = 1
  K = makeK(m, seed=8)
  y = simUVM(n, X, Z, K, seed=9)
  ans = checkUVM(y, X, Z, K)
  expect_true(all(ans$u[(n+1):m] != 0))
})

test_that("REJECT a kernel that is not positive semidefinite", {
  # Both decompositions check, so one record per level and two are tried
  set.seed(15)
  n = 10
  X = matrix(1, n, 1)
  y = matrix(rnorm(n), ncol=1)
  K = diag(n)
  K[1,1] = -1
  expect_error(solveUVM(y, X, diag(n), K), "positive semidefinite")
  K = diag(n/2)
  K[1,1] = -1
  expect_error(solveUVM(y, X, makeZ(n, n/2), K), "positive semidefinite")
})

test_that("REJECT inputs whose dimensions do not agree", {
  set.seed(16)
  n = 10
  y = matrix(rnorm(n), ncol=1)
  expect_error(solveUVM(y, matrix(1, n, 1), diag(n), diag(n+1)),
               "dimensions")
})

test_that("REJECT models that cannot be estimated", {
  set.seed(17)
  n = 12
  y = matrix(rnorm(n), ncol=1)
  X = matrix(1, n, 1)
  Z = makeZ(n, 4)
  # No residual degrees of freedom
  expect_error(solveUVM(y, matrix(rnorm(n*n), n, n), diag(n), diag(n)),
               "more records than fixed effects")
  # A random effect with no variance in any direction
  expect_error(solveUVM(y, X, Z, matrix(0, 4, 4)), "no positive eigenvalues")
  expect_error(solveUVM(y, X, matrix(0, n, 4), diag(4)), "is zero")
  # No fixed effect at all
  expect_error(solveUVM(y, matrix(0, n, 1), Z, diag(4)), "nonzero column")
  # An indefinite K whose negative eigenvalue is on a level with no record,
  # which Z*K*Z' cannot see
  expect_error(solveUVM(matrix(rnorm(3), ncol=1), matrix(1, 3, 1),
                        cbind(diag(3), 0), diag(c(1, 2, 3, -1))),
               "not positive semidefinite")
})

test_that("DROP dependent columns of X and keep the intercept", {
  # The intercept is placed last, after two indicators that add up to it, so
  # the greedy choice would keep both indicators if it did not consider the
  # intercept first. Dropping a dependent column must leave the fit on the
  # column space unchanged.
  set.seed(18)
  n = 30
  a = rep(c(1, 0), each=n/2)
  Z = makeZ(n, 6)
  K = crossprod(matrix(rnorm(36), 6, 6))/6
  y = matrix(2 + a + Z%*%rnorm(6) + rnorm(n), ncol=1)
  ref = solveUVM(y, cbind(1, a), Z, K)
  dep = solveUVM(y, cbind(a, 1-a, 1), Z, K)
  expect_equal(dep$Vu, ref$Vu, tolerance=1e-6)
  expect_equal(dep$Ve, ref$Ve, tolerance=1e-6)
  expect_equal(dep$LL, ref$LL, tolerance=1e-8)
  expect_equal(c(dep$u), c(ref$u), tolerance=1e-6)
  expect_equal(c(dep$beta), c(ref$beta[2], 0, ref$beta[1]), tolerance=1e-6)
  # A duplicated column
  dup = solveUVM(y, cbind(1, a, a), Z, K)
  expect_equal(c(dup$beta), c(ref$beta, 0), tolerance=1e-6)
  expect_equal(c(dup$u), c(ref$u), tolerance=1e-6)
})

test_that("SAME fit for an incidence Z as for an equivalent general Z", {
  # 30 records on 60 levels, half of them without a record, so Z*K*Z' is
  # formed on the records. For an incidence Z it is read out of K; for Z*Q
  # with K' = Q'*K*Q, the same model, it comes from the dense products. The
  # BLUPs of the second are Q' times those of the first.
  set.seed(20)
  n = 30
  m = 60
  lev = c(1:20, sample(20, n-20, replace=TRUE))
  Z = matrix(0, n, m)
  Z[cbind(1:n, lev)] = 1
  B = matrix(rnorm(m*80), m)
  K = tcrossprod(B)/80
  y = matrix(1 + Z%*%(t(chol(K))%*%rnorm(m)) + rnorm(n), ncol=1)
  X = matrix(1, n, 1)
  Q = qr.Q(qr(matrix(rnorm(m*m), m)))
  a = solveUVM(y, X, Z, K)
  b = solveUVM(y, X, Z%*%Q, t(Q)%*%K%*%Q)
  expect_equal(a$Vu, b$Vu, tolerance=1e-6)
  expect_equal(a$Ve, b$Ve, tolerance=1e-6)
  expect_equal(a$LL, b$LL, tolerance=1e-8)
  expect_equal(c(a$beta), c(b$beta), tolerance=1e-6)
  expect_equal(c(a$u), c(Q%*%b$u), tolerance=1e-6)
  # Levels without a record still get BLUPs through K
  expect_true(all(a$u[41:60]!=0))
})

test_that("SHIFT a large mean or covariate offset leaves the fit unchanged", {
  # A response with a mean of 1e8 would lose its variation to rounding if
  # the residual sum of squares were found by subtracting sums of squares,
  # and a covariate of 1e6+x would look like a copy of the intercept if
  # independence were judged on X'X. Neither may change the variances, and
  # the fixed effects must move exactly as the reparameterization implies.
  set.seed(21)
  n = 40
  Z = makeZ(n, 8)
  K = crossprod(matrix(rnorm(64), 8, 8))/8
  x = rnorm(n)
  y = matrix(1 + 0.5*x + Z%*%rnorm(8) + rnorm(n), ncol=1)
  ref = solveUVM(y, cbind(1, x), Z, K)
  big = solveUVM(y + 1e8, cbind(1, x), Z, K)
  expect_equal(big$Vu, ref$Vu, tolerance=1e-6)
  expect_equal(big$Ve, ref$Ve, tolerance=1e-6)
  expect_equal(c(big$beta) - c(1e8, 0), c(ref$beta), tolerance=1e-6,
               scale=1)
  off = solveUVM(y, cbind(1, 1e6 + x), Z, K)
  expect_equal(off$Vu, ref$Vu, tolerance=1e-6)
  expect_equal(off$Ve, ref$Ve, tolerance=1e-6)
  expect_equal(off$beta[2], ref$beta[2], tolerance=1e-6)
  # The same for solveRRBLUP, on both of its decompositions
  for(m in c(10, 60)){
    M = matrix(sample(0:2, n*m, replace=TRUE), n, m)
    yM = matrix(1 + 0.5*x + M%*%rnorm(m, sd=0.2) + rnorm(n), ncol=1)
    a = solveRRBLUP(yM, cbind(1, x), M)
    b = solveRRBLUP(yM + 1e8, cbind(1, 1e6 + x), M)
    expect_equal(b$Vu, a$Vu, tolerance=1e-6)
    expect_equal(b$Ve, a$Ve, tolerance=1e-6)
    expect_equal(b$beta[2], a$beta[2], tolerance=1e-6)
    expect_equal(c(b$u), c(a$u), tolerance=1e-6)
  }
})

# 3. solveRRBLUP ----------------------------------------------------------

test_that("AGREE solveRRBLUP with solveUVM on a linear kernel", {
  # RR-BLUP is the model Z = M, K = I, so the variances must match and the
  # marker effects are the BLUPs
  set.seed(11)
  n = 25
  m = 12
  X = cbind(1, rnorm(n))
  M = matrix(sample(0:2, n*m, replace=TRUE), n, m)
  y = matrix(M%*%rnorm(m, sd=0.3) + rnorm(n), ncol=1)
  a = solveRRBLUP(y, X, M)
  b = solveUVM(y, X, M, diag(m))
  expect_equal(a$Vu, b$Vu, tolerance=1e-6)
  expect_equal(a$Ve, b$Ve, tolerance=1e-6)
  expect_equal(c(a$beta), c(b$beta), tolerance=1e-6)
  expect_equal(c(a$u), c(b$u), tolerance=1e-6)
})

test_that("DROP dependent columns of X in solveRRBLUP and keep the intercept", {
  # Both decompositions, fewer markers than records and more
  set.seed(19)
  n = 30
  a = rep(c(1, 0), each=n/2)
  for(m in c(10, 50)){
    M = matrix(sample(0:2, n*m, replace=TRUE), n, m)
    y = matrix(1 + a + M%*%rnorm(m, sd=0.2) + rnorm(n), ncol=1)
    ref = solveRRBLUP(y, cbind(1, a), M)
    dep = solveRRBLUP(y, cbind(a, 1-a, 1), M)
    expect_equal(dep$Vu, ref$Vu, tolerance=1e-6)
    expect_equal(dep$Ve, ref$Ve, tolerance=1e-6)
    expect_equal(c(dep$u), c(ref$u), tolerance=1e-6)
    expect_equal(c(dep$beta), c(ref$beta[2], 0, ref$beta[1]), tolerance=1e-6)
  }
  expect_error(solveRRBLUP(y, matrix(0, n, 2), M), "nonzero column")
  expect_error(solveRRBLUP(y, matrix(rnorm(n*n), n, n), M),
               "rank of X must be smaller")
})

test_that("RUN RRBLUP end to end", {
  set.seed(12)
  founderPop = quickHaplo(nInd=40, nChr=2, segSites=20)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=5)
  SP$addSnpChip(nSnpPerChr=10)
  SP$setVarE(h2=0.5)
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)
  ans = RRBLUP(pop, simParam=SP)
  expect_true(all(is.finite(ans@gv[[1]]@addEff)))
  expect_gt(ans@Vu[1,1], 0)
  expect_gt(ans@Ve[1,1], 0)
  pop = setEBV(pop, ans, simParam=SP)
  expect_true(all(is.finite(pop@ebv)))
})
