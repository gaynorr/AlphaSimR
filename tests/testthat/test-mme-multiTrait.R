context("mme-multiTrait")

# The multiple trait solvers are checked against what must hold for any
# correct REML fit rather than against stored numbers. A dense reference
# builds the full covariance of vec(Y) and maximizes the restricted
# log-likelihood numerically, so the solvers have to reach at least its
# maximum. At whatever covariances they return, the fixed and random effects
# have to solve the mixed model equations.

# ---------------------------------------------------------------------------
# 1. Dense reference
# ---------------------------------------------------------------------------

# Restricted log-likelihood of vec(Y) with the convention the solvers use,
# which leaves out the log|X'X| term
denseLL = function(Vu, Ve, Y, X, Kn){
  n = nrow(Y)
  nt = ncol(Y)
  y = c(Y)
  Xb = kronecker(diag(nt), X)
  V = kronecker(Vu, Kn) + kronecker(Ve, diag(n))
  cV = chol(V)
  Vi = chol2inv(cV)
  XVX = crossprod(Xb, Vi %*% Xb)
  cX = chol(XVX)
  ViX = Vi %*% Xb
  P = Vi - ViX %*% chol2inv(cX) %*% t(ViX)
  p = ncol(Xb)
  -0.5*((n*nt-p)*log(2*pi) + 2*sum(log(diag(cV))) +
          2*sum(log(diag(cX))) + drop(crossprod(y, P %*% y)))
}

# Builds a covariance matrix from the elements of its lower Cholesky factor
lowerToCov = function(par, nt){
  L = matrix(0, nt, nt)
  L[lower.tri(L, diag=TRUE)] = par
  tcrossprod(L)
}

# Maximizes denseLL over Cholesky factors of Vu and Ve from a start that owes
# nothing to the solvers
denseMax = function(Y, X, Kn){
  nt = ncol(Y)
  nPar = nt*(nt+1)/2
  v = apply(Y, 2, var)/2
  L0 = t(chol(diag(v, nt)))
  start = c(L0[lower.tri(L0, diag=TRUE)], L0[lower.tri(L0, diag=TRUE)])
  obj = function(par){
    Vu = lowerToCov(par[1:nPar], nt)
    Ve = lowerToCov(par[nPar+(1:nPar)], nt)
    ans = tryCatch(denseLL(Vu, Ve, Y, X, Kn), error=function(e) -1e10)
    -ans
  }
  best = optim(start, obj, method="BFGS",
               control=list(maxit=5000, reltol=1e-12))
  best = optim(best$par, obj, method="Nelder-Mead",
               control=list(maxit=20000, reltol=1e-14))
  best = optim(best$par, obj, method="BFGS",
               control=list(maxit=5000, reltol=1e-14))
  -best$value
}

# ---------------------------------------------------------------------------
# 2. Simulated data
# ---------------------------------------------------------------------------

# Two traits with correlated genetic and error effects on a kernel Kn
simTraits = function(Kn, X, seed, nt=2){
  set.seed(seed)
  n = nrow(Kn)
  Vu = matrix(c(1, 0.5, 0.5, 1), 2)[1:nt, 1:nt, drop=FALSE]
  Ve = matrix(c(1, -0.3, -0.3, 1.5), 2)[1:nt, 1:nt, drop=FALSE]
  e = eigen(Kn, symmetric=TRUE)
  Kh = e$vectors %*% diag(sqrt(pmax(e$values, 0))) %*% t(e$vectors)
  G = Kh %*% matrix(rnorm(n*nt), n, nt) %*% chol(Vu)
  E = matrix(rnorm(n*nt), n, nt) %*% chol(Ve)
  B = matrix(seq_len(ncol(X)*nt), ncol(X), nt)
  X %*% B + G + E
}

markers = function(n, m, seed){
  set.seed(seed)
  M = matrix(rbinom(n*m, 2, 0.5), n, m)
  scale(M, scale=FALSE)
}

# Henderson's mixed model equations for vec(Y). They need the inverses of K,
# Vu and Ve, so they are checked only where the estimates are well inside the
# parameter space. With few records the REML estimate of a correlation is
# often at or near one, which leaves Vu or Ve close to singular.
expectMME = function(ans, Y, X, Z, K){
  n = nrow(Y)
  nt = ncol(Y)
  y = c(Y)
  Xb = kronecker(diag(nt), X)
  Zb = kronecker(diag(nt), Z)
  Ri = kronecker(solve(ans$Ve), diag(n))
  Gi = kronecker(solve(ans$Vu), solve(K))
  C = rbind(cbind(crossprod(Xb, Ri %*% Xb), crossprod(Xb, Ri %*% Zb)),
            cbind(crossprod(Zb, Ri %*% Xb), crossprod(Zb, Ri %*% Zb) + Gi))
  rhs = c(crossprod(Xb, Ri %*% y), crossprod(Zb, Ri %*% y))
  expect_equal(c(C %*% c(ans$beta, ans$u)), rhs, tolerance=1e-6)
}

# The same equations written on the records, which hold whenever V is
# positive definite, including when K, Vu or Ve is singular:
# beta is the generalized least squares solution and u = (Vu (x) K*Z')*V^-1*r
expectGLS = function(ans, Y, X, Z, K){
  n = nrow(Y)
  nt = ncol(Y)
  y = c(Y)
  Xb = kronecker(diag(nt), X)
  V = kronecker(ans$Vu, Z %*% K %*% t(Z)) + kronecker(ans$Ve, diag(n))
  Vi = solve(V)
  beta = solve(crossprod(Xb, Vi %*% Xb), crossprod(Xb, Vi %*% y))
  expect_equal(c(ans$beta), c(beta), tolerance=1e-6)
  u = kronecker(ans$Vu, K %*% t(Z)) %*% Vi %*% (y - Xb %*% beta)
  expect_equal(c(ans$u), c(u), tolerance=1e-6)
}

# Checks shared by every fit
expectFit = function(ans, Y, X, Kn, mLevels){
  nt = ncol(Y)
  expect_equal(dim(ans$Vu), c(nt, nt))
  expect_equal(dim(ans$Ve), c(nt, nt))
  expect_equal(dim(ans$beta), c(ncol(X), nt))
  expect_equal(dim(ans$u), c(mLevels, nt))
  expect_equal(ans$Vu, t(ans$Vu))
  expect_equal(ans$Ve, t(ans$Ve))
  expect_true(all(eigen(ans$Ve, symmetric=TRUE)$values > 0))
  expect_true(all(eigen(ans$Vu, symmetric=TRUE)$values > -1e-8))
  # The reported likelihood is the dense one at the returned covariances
  expect_equal(ans$LL, denseLL(ans$Vu, ans$Ve, Y, X, Kn), tolerance=1e-8)
  # and is no lower than the numerical maximum
  expect_gte(ans$LL, denseMax(Y, X, Kn) - 1e-4)
}

# ---------------------------------------------------------------------------
# 3. solveMVM
# ---------------------------------------------------------------------------

test_that("MVM a full rank kernel reaches the REML maximum", {
  n = 40
  M = markers(n, 60, 9101)
  K = tcrossprod(M)/60 + diag(0.05, n)
  X = matrix(1, n, 1)
  Y = simTraits(K, X, 9102)
  Z = diag(n)
  ans = solveMVM(Y, X, Z, K)
  expectFit(ans, Y, X, K, n)
  expectGLS(ans, Y, X, Z, K)
  expect_true(ans$iter >= 1)
})

test_that("MVM repeated records per level reach the REML maximum", {
  # Fifteen levels with three records each, so Z*K*Z' has rank fifteen
  m = 15
  M = markers(m, 40, 9201)
  K = tcrossprod(M)/40 + diag(0.05, m)
  Z = kronecker(diag(m), matrix(1, 3, 1))
  n = nrow(Z)
  X = cbind(1, rep(c(-1, 0, 1), m))
  Kn = Z %*% K %*% t(Z)
  Y = simTraits(Kn, X, 9202)
  ans = solveMVM(Y, X, Z, K)
  expectFit(ans, Y, X, Kn, m)
  expectGLS(ans, Y, X, Z, K)
})

test_that("MVM a singular K is handled without breakdown", {
  # K has rank eight on twenty levels, and there are more records than levels
  m = 20
  M = markers(m, 8, 9301)
  K = tcrossprod(M)/8
  Z = kronecker(diag(m), matrix(1, 2, 1))
  n = nrow(Z)
  X = matrix(1, n, 1)
  Kn = Z %*% K %*% t(Z)
  Y = simTraits(Kn, X, 9302)
  ans = solveMVM(Y, X, Z, K)
  expect_true(all(is.finite(unlist(ans))))
  expectFit(ans, Y, X, Kn, m)
  expectGLS(ans, Y, X, Z, K)
})

test_that("MVM a singular K is handled when there are more levels than records", {
  m = 30
  M = markers(m, 6, 9351)
  K = tcrossprod(M)/6
  Z = diag(m)[sample(m, 25), ]
  X = matrix(1, 25, 1)
  Kn = Z %*% K %*% t(Z)
  Y = simTraits(Kn, X, 9352)
  ans = solveMVM(Y, X, Z, K)
  expect_true(all(is.finite(unlist(ans))))
  expectFit(ans, Y, X, Kn, m)
  expectGLS(ans, Y, X, Z, K)
})

test_that("MVM estimates inside the parameter space solve Henderson's equations", {
  m = 50
  M = markers(m, 40, 9251)
  K = tcrossprod(M)/40 + diag(0.05, m)
  Z = kronecker(diag(m), matrix(1, 3, 1))
  X = cbind(1, rep(c(-1, 0, 1), m))
  Kn = Z %*% K %*% t(Z)
  Y = simTraits(Kn, X, 9252)
  ans = solveMVM(Y, X, Z, K)
  expect_lt(kappa(ans$Vu, exact=TRUE), 100)
  expect_lt(kappa(ans$Ve, exact=TRUE), 100)
  expect_equal(ans$LL, denseLL(ans$Vu, ans$Ve, Y, X, Kn), tolerance=1e-8)
  expectMME(ans, Y, X, Z, K)
  expectGLS(ans, Y, X, Z, K)
})

test_that("MVM an incidence Z gives the same fit as an equivalent general Z", {
  # 30 records on 60 levels, so Z*K*Z' is formed on the records. For an
  # incidence Z it is read out of K; for Z*Q with K' = Q'*K*Q, the same
  # model, it comes from the dense products.
  set.seed(9941)
  n = 30
  m = 60
  lev = c(1:20, sample(20, n-20, replace=TRUE))
  Z = matrix(0, n, m)
  Z[cbind(1:n, lev)] = 1
  B = matrix(rnorm(m*80), m)
  K = tcrossprod(B)/80
  Y = simTraits(Z%*%K%*%t(Z), matrix(1, n, 1), 9942)
  X = matrix(1, n, 1)
  Q = qr.Q(qr(matrix(rnorm(m*m), m)))
  a = solveMVM(Y, X, Z, K)
  b = solveMVM(Y, X, Z%*%Q, t(Q)%*%K%*%Q)
  expect_equal(a$Vu, b$Vu, tolerance=1e-6)
  expect_equal(a$Ve, b$Ve, tolerance=1e-6)
  expect_equal(a$LL, b$LL, tolerance=1e-8)
  expect_equal(a$beta, b$beta, tolerance=1e-6)
  expect_equal(a$u, Q%*%b$u, tolerance=1e-6)
})

test_that("MVM an indefinite K is rejected on both decompositions", {
  # Clipping its negative eigenvalue would fit one kernel while the BLUPs
  # used another
  set.seed(9943)
  K = diag(c(1, 2, 3, -0.1))
  X = matrix(1, 4, 1)
  Y = matrix(rnorm(8), 4, 2)
  # As many levels as records: the kernel on the records is decomposed
  expect_error(solveMVM(Y, X, diag(4), K), "positive semidefinite")
  # Fewer levels than records: K itself is factored
  n = 12
  Z = matrix(0, n, 4)
  Z[cbind(1:n, rep(1:4, 3))] = 1
  expect_error(solveMVM(matrix(rnorm(2*n), n, 2), matrix(1, n, 1), Z, K),
               "positive semidefinite")
  # A negative eigenvalue on a level with no record, which Z*K*Z' cannot see
  expect_error(solveMVM(matrix(rnorm(6), 3, 2), matrix(1, 3, 1),
                        cbind(diag(3), 0), diag(c(1, 2, 3, -1))),
               "positive semidefinite")
})

test_that("RRBLUPMV marker effects are recovered when the residual variance is near zero", {
  # Nearly noiseless traits put Ve close to zero. The part of P*y outside the
  # markers' span contributes nothing to the effects, and once added rounding
  # amplified by Ve^-1 that left them wrong by several percent.
  n = 300
  m = 40
  M = markers(n, m, 9946)
  set.seed(9947)
  A = matrix(rnorm(2*m, sd=0.2), m, 2)
  Y = 1 + M %*% A + matrix(rnorm(2*n, sd=1e-6), n, 2)
  ans = solveRRBLUPMV(Y, matrix(1, n, 1), M)
  expect_lt(max(abs(ans$u - A))/max(abs(A)), 1e-3)
})

test_that("MVM the BLUPs use the kernel that was fitted", {
  # The third level's eigenvalue, -1e-9, is taken as rounding and dropped
  # from the fitted kernel, so the BLUPs must give that level no effect
  set.seed(9948)
  n = 30
  Z = matrix(0, n, 3)
  Z[cbind(1:n, rep(1:3, 10))] = 1
  K = diag(c(1, 2, -1e-9))
  Y = matrix(rnorm(2*n), n, 2) + Z %*% matrix(rnorm(6), 3, 2)
  ans = solveMVM(Y, matrix(1, n, 1), Z, K)
  expect_lt(max(abs(ans$u[3, ])), 1e-12)
})

test_that("MVM a large mean does not change the covariance matrices", {
  n = 30
  m = 10
  M = markers(n, m, 9944)
  X = matrix(1, n, 1)
  Y = simTraits(tcrossprod(M)/m, X, 9945)
  a = solveRRBLUPMV(Y, X, M)
  b = solveRRBLUPMV(Y + 1e8, X, M)
  expect_equal(b$Vu, a$Vu, tolerance=1e-6)
  expect_equal(b$Ve, a$Ve, tolerance=1e-6)
  expect_equal(c(b$beta) - 1e8, c(a$beta), tolerance=1e-6, scale=1)
})

# ---------------------------------------------------------------------------
# 4. solveRRBLUPMV
# ---------------------------------------------------------------------------

test_that("RRBLUPMV fewer markers than records reach the REML maximum", {
  n = 40
  m = 10
  M = markers(n, m, 9401)
  X = matrix(1, n, 1)
  Kn = tcrossprod(M)
  Y = simTraits(Kn/m, X, 9402)
  ans = solveRRBLUPMV(Y, X, M)
  expectFit(ans, Y, X, Kn, m)
  expectGLS(ans, Y, X, M, diag(m))
})

test_that("RRBLUPMV more markers than records reach the REML maximum", {
  n = 30
  m = 50
  M = markers(n, m, 9501)
  X = cbind(1, rnorm(n))
  Kn = tcrossprod(M)
  Y = simTraits(Kn/m, X, 9502)
  ans = solveRRBLUPMV(Y, X, M)
  expectFit(ans, Y, X, Kn, m)
  expectGLS(ans, Y, X, M, diag(m))
})

test_that("RRBLUPMV estimates inside the parameter space solve Henderson's equations", {
  n = 150
  m = 30
  M = markers(n, m, 101)
  X = cbind(1, rnorm(n))
  Kn = tcrossprod(M)
  Y = simTraits(Kn/m, X, 201)
  ans = solveRRBLUPMV(Y, X, M)
  expect_lt(kappa(ans$Vu, exact=TRUE), 100)
  expect_lt(kappa(ans$Ve, exact=TRUE), 100)
  expect_equal(ans$LL, denseLL(ans$Vu, ans$Ve, Y, X, Kn), tolerance=1e-8)
  expectMME(ans, Y, X, M, diag(m))
  expectGLS(ans, Y, X, M, diag(m))
})

test_that("RRBLUPMV and MVM agree when K is the identity on the markers", {
  n = 35
  m = 12
  M = markers(n, m, 9601)
  X = matrix(1, n, 1)
  Y = simTraits(tcrossprod(M)/m, X, 9602)
  a = solveRRBLUPMV(Y, X, M)
  b = solveMVM(Y, X, M, diag(m))
  expect_equal(a$Vu, b$Vu, tolerance=1e-5)
  expect_equal(a$Ve, b$Ve, tolerance=1e-5)
  expect_equal(a$u, b$u, tolerance=1e-5)
  expect_equal(a$LL, b$LL, tolerance=1e-8)
})

test_that("RRBLUPMV a single trait agrees with a univariate REML search", {
  n = 40
  m = 25
  M = markers(n, m, 9701)
  X = matrix(1, n, 1)
  Kn = tcrossprod(M)
  Y = simTraits(Kn/m, X, 9702, nt=1)
  ans = solveRRBLUPMV(Y, X, M)
  # The univariate likelihood has a closed form once log(Vu) and log(Ve) are
  # given, so a two parameter search is an exact reference
  obj = function(par){
    -denseLL(matrix(exp(par[1])), matrix(exp(par[2])), Y, X, Kn)
  }
  best = optim(c(0, 0), obj, method="BFGS",
               control=list(maxit=5000, reltol=1e-14))
  expect_gte(ans$LL, -best$value - 1e-6)
  expect_equal(c(ans$Vu), exp(best$par[1]), tolerance=1e-3)
  expect_equal(c(ans$Ve), exp(best$par[2]), tolerance=1e-3)
  expectMME(ans, Y, X, M, diag(m))
})

test_that("RRBLUPMV the fit does not depend on the order of the traits", {
  n = 40
  m = 15
  M = markers(n, m, 9801)
  X = matrix(1, n, 1)
  Y = simTraits(tcrossprod(M)/m, X, 9802)
  a = solveRRBLUPMV(Y, X, M)
  b = solveRRBLUPMV(Y[, 2:1], X, M)
  expect_equal(a$Vu, b$Vu[2:1, 2:1], tolerance=1e-5)
  expect_equal(a$Ve, b$Ve[2:1, 2:1], tolerance=1e-5)
  expect_equal(a$LL, b$LL, tolerance=1e-8)
})

test_that("RRBLUPMV warns when maxIter is reached", {
  n = 30
  m = 10
  M = markers(n, m, 9901)
  X = matrix(1, n, 1)
  Y = simTraits(tcrossprod(M)/m, X, 9902)
  expect_warning(ans <- solveRRBLUPMV(Y, X, M, maxIter=2L), "converge")
  expect_equal(ans$iter, 2L)
})

test_that("INPUT malformed arguments are refused with a message", {
  n = 20
  m = 8
  M = markers(n, m, 9911)
  X = matrix(1, n, 1)
  Y = simTraits(tcrossprod(M)/m, X, 9912)
  expect_error(solveRRBLUPMV(Y, X[-1, , drop=FALSE], M), "same number of rows")
  expect_error(solveRRBLUPMV(Y[, 0, drop=FALSE], X, M), "at least one column")
  expect_error(solveRRBLUPMV(Y, X[, 0, drop=FALSE], M), "at least one column")
  expect_error(solveRRBLUPMV(Y, matrix(rnorm(n*n), n, n), M),
               "rank of X must be smaller")
  expect_error(solveRRBLUPMV(Y, matrix(0, n, 2), M), "nonzero column")
  Yna = Y
  Yna[3, 1] = NA
  expect_error(solveRRBLUPMV(Yna, X, M), "missing")
  expect_error(solveRRBLUPMV(Y, X, M[-1, ]), "Y and M")
  expect_error(solveMVM(Y, X, diag(n)[-1, ], diag(n)), "Y and Z")
  expect_error(solveMVM(Y, X, diag(n), diag(n - 1)), "K must be square")
})

test_that("INPUT a kernel with no positive eigenvalue is refused", {
  n = 20
  X = matrix(1, n, 1)
  Y = matrix(rnorm(2*n), n, 2)
  # Both ways of rotating the records: fewer markers than records, and more
  expect_error(solveRRBLUPMV(Y, X, matrix(0, n, 5)), "no positive eigenvalue")
  expect_error(solveRRBLUPMV(Y, X, matrix(0, n, 30)), "no positive eigenvalue")
})

test_that("DROP dependent columns of X and keep the intercept", {
  # The intercept is placed last, after two indicators that add up to it, so
  # the greedy choice would keep both indicators if it did not consider the
  # intercept first. Dropping a dependent column must leave the fit on the
  # column space unchanged. Both ways of rotating the records are used.
  n = 30
  a = rep(c(1, 0), each=n/2)
  for(m in c(10, 50)){
    M = markers(n, m, 9921+m)
    Y = simTraits(tcrossprod(M)/m, cbind(1, a), 9922+m)
    ref = solveRRBLUPMV(Y, cbind(1, a), M)
    dep = solveRRBLUPMV(Y, cbind(a, 1-a, 1), M)
    expect_equal(dep$Vu, ref$Vu, tolerance=1e-6)
    expect_equal(dep$Ve, ref$Ve, tolerance=1e-6)
    expect_equal(dep$LL, ref$LL, tolerance=1e-8)
    expect_equal(dep$u, ref$u, tolerance=1e-6)
    expect_equal(dep$beta, rbind(ref$beta[2, ], 0, ref$beta[1, ]),
                 tolerance=1e-6)
    mvm = solveMVM(Y, cbind(a, 1-a, 1), diag(n), tcrossprod(M)/m)
    expect_equal(mvm$beta[2, ], c(0, 0))
    expect_equal(mvm$LL, solveMVM(Y, cbind(1, a), diag(n),
                                  tcrossprod(M)/m)$LL, tolerance=1e-8)
  }
})

# ---------------------------------------------------------------------------
# 5. End to end
# ---------------------------------------------------------------------------

test_that("RRBLUP on two traits routes through the multiple trait solver", {
  skip_on_cran()
  set.seed(9951)
  founderPop = quickHaplo(nInd=60, nChr=2, segSites=30)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$restrSegSites(minQtlPerChr=10, minSnpPerChr=15, overlap=FALSE)
  SP$addTraitA(nQtlPerChr=10, mean=c(0,0), var=c(1,1),
               corA=matrix(c(1,0.5,0.5,1), nrow=2))
  SP$addSnpChip(nSnpPerChr=15)
  SP$setVarE(h2=c(0.5,0.5))
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)

  ans = RRBLUP(pop, traits=1:2, simParam=SP)

  # RRBLUP codes a diploid's genotypes as -1, 0 and 1 and has a single fixed
  # effect level here, so the same fit can be made directly
  M = pullSnpGeno(pop, simParam=SP) - 1
  ref = solveRRBLUPMV(pheno(pop), matrix(1, nInd(pop), 1), M)
  expect_equal(ans@Vu, ref$Vu)
  expect_equal(ans@Ve, ref$Ve)
  for(i in 1:2){
    expect_equal(c(ans@bv[[i]]@addEff), c(ref$u[, i]))
    expect_equal(ans@gv[[i]]@intercept, ref$beta[1, i])
    expect_equal(ans@bv[[i]]@intercept, -sum(colMeans(M)*ref$u[, i]))
  }

  pred = setEBV(pop, ans, simParam=SP)
  expect_equal(ncol(ebv(pred)), 2L)
  expect_true(all(is.finite(c(ebv(pred)))))
})
