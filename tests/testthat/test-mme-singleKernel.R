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
  # Fixed effects that are not identifiable
  expect_error(solveUVM(y, cbind(1, 1, rnorm(n)), Z, diag(4)), "full column rank")
})

# 3. solveRRBLUP ----------------------------------------------------------

test_that("AGREE solveRRBLUP with solveRRBLUP2 on either side", {
  set.seed(10)
  for(m in c(8, 40)){
    n = 25
    X = cbind(1, rnorm(n))
    M = matrix(sample(0:2, n*m, replace=TRUE), n, m)
    y = matrix(M%*%rnorm(m, sd=0.3) + rnorm(n), ncol=1)
    expect_identical(solveRRBLUP(y, X, M), solveRRBLUP2(y, X, M))
  }
})

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
