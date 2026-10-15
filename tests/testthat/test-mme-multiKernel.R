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

test_that("MAXITER reaching the iteration limit gives a warning", {
  d = simTwoKernels(seed=9110)
  expect_warning(solveRRBLUPMK(d$y, d$X, list(d$M1, d$M2), maxIter=1L),
                 "did not converge")
  expect_warning(solveMKM(d$y, d$X, list(diag(d$n)),
                          list(tcrossprod(d$M1)), maxIter=1L),
                 "did not converge")
})

test_that("INPUT mismatched dimensions are rejected", {
  d = simTwoKernels(seed=9111)
  expect_error(solveRRBLUPMK(d$y, d$X, list(d$M1[-1, ])), "one row per record")
  expect_error(solveMKM(d$y, d$X, list(diag(d$n)), list(diag(3))),
               "square")
  expect_error(solveMKM(d$y, d$X, list(diag(d$n)), list()), "same length")
})

# ---------------------------------------------------------------------------
# 5. The genomic selection models that call solveRRBLUPMK
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
