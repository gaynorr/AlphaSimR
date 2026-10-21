context("pedBLUP")

# The internal fit is checked against two references that need no sparse
# algebra. The BLUP is written in its V^-1 form, u = Vu*A*Z'*V^-1*(y-X*b),
# which holds even when A is singular, as it is with doubled haploids of
# inbred parents. REML is checked against solveUVM's dense fit on the
# relationships among the phenotyped individuals, which is the same
# likelihood because ancestors without records integrate out.

# 1. Helpers ----

denseBLUP = function(y, X, Z, A, Vu, Ve){
  V = Vu*Z %*% A %*% t(Z) + Ve*diag(length(y))
  Vi = solve(V)
  b = solve(t(X) %*% Vi %*% X, t(X) %*% Vi %*% y)
  u = Vu*A %*% t(Z) %*% Vi %*% (y - X %*% b)
  list(b=c(b), u=c(u))
}

# A simulation with crosses, doubled haploids, selection and four levels of
# a fixed effect
blupSim = function(seed, ploidy=2L, nDH=10){
  set.seed(seed)
  founderPop = quickHaplo(nInd=40, nChr=2, segSites=50, ploidy=ploidy)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackPed(TRUE)
  SP$setSexes("no")
  SP$addTraitA(nQtlPerChr=20)
  SP$setVarE(h2=0.4)
  pop = newPop(founderPop, simParam=SP)
  train = pop
  for(g in 1:4){
    pop = randCross(selectInd(pop, 20, simParam=SP), nCrosses=60,
                    simParam=SP)
    if(nDH>0){
      pop = c(pop, makeDH(pop[seq_len(nDH)], simParam=SP))
    }
    pop = setPheno(pop, fixEff=as.integer(g), simParam=SP)
    train = c(train, pop)
  }
  list(SP=SP, train=train)
}

# One column per level of the fixed effects, which model.matrix cannot
# give for a single level
fixedDesign = function(fixEff){
  f = factor(fixEff)
  outer(as.integer(f), seq_len(nlevels(f)), "==")*1
}

# The dense design of a fit: X, Z over the pruned pedigree, and its A
denseDesign = function(sim, train=sim$train){
  ped = AlphaSimR:::.prunePed(train@iid, train@ploidy, sim$SP)
  X = fixedDesign(train@fixEff)
  Z = matrix(0, train@nInd, length(ped$iid))
  Z[cbind(seq_len(train@nInd), ped$row)] = 1
  A = AlphaSimR:::calcPedA(ped$mother, ped$father, ped$isDH,
                           ploidy=train@ploidy)
  list(ped=ped, X=X, Z=Z, A=A)
}

fit = function(sim, train=sim$train, ...){
  AlphaSimR:::.fitPedBLUP(train@pheno[,1], train@iid, train@fixEff,
                          train@ploidy, simParam=sim$SP, ...)
}

# 2. Given variance components ----

test_that("BLUP matches Henderson's BLUP with given variance components", {
  sim = blupSim(9101)
  d = denseDesign(sim)
  ans = fit(sim, Vu=0.7, Ve=1.3)
  ref = denseBLUP(sim$train@pheno[,1], d$X, d$Z, d$A, 0.7, 1.3)
  expect_identical(ans$iid, d$ped$iid)
  expect_equal(ans$ebv, ref$u, tolerance=1e-8)
  expect_equal(unname(ans$fixEff), ref$b, tolerance=1e-8)
  expect_equal(names(ans$fixEff), as.character(1:4))
  expect_equal(c(ans$Vu, ans$Ve), c(0.7, 1.3))
})

test_that("BLUP handles clones and F1s of inbred parents", {
  # Inbred founders make doubled haploids of F1s and F1s of doubled
  # haploids clones or exact averages of their parents
  set.seed(9102)
  founderPop = quickHaplo(nInd=20, nChr=1, segSites=50, inbred=TRUE)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackPed(TRUE)
  SP$setSexes("no")
  SP$addTraitA(nQtlPerChr=20)
  SP$setVarE(h2=0.5)
  pop = newPop(founderPop, simParam=SP)
  F1 = randCross(pop, nCrosses=30, simParam=SP)
  DH = makeDH(F1, nDH=2, simParam=SP)
  F1b = randCross(DH, nCrosses=30, simParam=SP)
  train = c(pop, F1, DH, F1b)
  sim = list(SP=SP, train=train)
  d = denseDesign(sim)
  ans = fit(sim, Vu=1, Ve=1)
  ref = denseBLUP(train@pheno[,1], d$X, d$Z, d$A, 1, 1)
  expect_equal(ans$ebv, ref$u, tolerance=1e-8)
})

test_that("REPEAT records of one individual share its breeding value", {
  sim = blupSim(9103)
  tr = sim$train
  rep2 = tr[sample.int(tr@nInd, 50)]
  rep2 = setPheno(rep2, simParam=sim$SP)
  train = c(tr, rep2)
  d = denseDesign(sim, train)
  ans = fit(sim, train, Vu=0.8, Ve=1.1)
  ref = denseBLUP(train@pheno[,1], d$X, d$Z, d$A, 0.8, 1.1)
  expect_equal(ans$ebv, ref$u, tolerance=1e-8)
})

test_that("POLYPLOID tetraploid BLUP matches Henderson's BLUP", {
  sim = blupSim(9104, ploidy=4L, nDH=0)
  d = denseDesign(sim)
  ans = fit(sim, Vu=0.6, Ve=1.2)
  ref = denseBLUP(sim$train@pheno[,1], d$X, d$Z, d$A, 0.6, 1.2)
  expect_equal(ans$ebv, ref$u, tolerance=1e-8)
})

# 3. REML ----

test_that("REML matches solveUVM's dense REML", {
  sim = blupSim(9105)
  tr = sim$train
  ans = fit(sim)
  K = pedigreeA(tr, simParam=sim$SP)
  X = fixedDesign(tr@fixEff)
  ref = solveUVM(tr@pheno, X, diag(tr@nInd), unname(K))
  # The search stops within 1e-3 of the best log(lambda), so the
  # variance components agree to about that
  expect_equal(ans$Vu, c(ref$Vu), tolerance=2e-3)
  expect_equal(ans$Ve, c(ref$Ve), tolerance=2e-3)
  expect_false(ans$boundary)
  # The BLUP at the estimates matches the dense BLUP at them
  d = denseDesign(sim)
  den = denseBLUP(tr@pheno[,1], d$X, d$Z, d$A, ans$Vu, ans$Ve)
  expect_equal(ans$ebv, den$u, tolerance=1e-8)
})

test_that("MULTIMODAL the REML search keeps the best it finds in its range", {
  # Two records of an unrelated doubled haploid and seven records from each
  # of two unrelated full-sib families. With the intercept projected out,
  # the relationships among the records have eigenvalues 0, 0.5 and 4, and
  # responses built from those eigenspaces give a likelihood with more
  # than one optimum. The search covers two orders of magnitude either
  # side of its starting value and is not meant to find a global optimum,
  # but it must stay in that range and do at least as well as the best of
  # its grid points.
  set.seed(9210)
  founderPop = quickHaplo(nInd=5, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackPed(TRUE)
  SP$setSexes("no")
  pop = newPop(founderPop, simParam=SP)
  dh = makeDH(pop[5], simParam=SP)
  famA = makeCross(pop, cbind(1L, 2L), nProgeny=7, simParam=SP)
  famB = makeCross(pop, cbind(3L, 4L), nProgeny=7, simParam=SP)
  iid = c(dh@iid, dh@iid, famA@iid, famB@iid)
  n = length(iid)
  ped = AlphaSimR:::.prunePed(iid, 2L, SP)
  A = AlphaSimR:::calcPedA(ped$mother, ped$father, ped$isDH)
  K = A[ped$row, ped$row]
  P = diag(n) - matrix(1/n, n, n)
  e = eigen(P %*% K %*% P, symmetric=TRUE)
  v4 = e$vectors[, which.min(abs(e$values-4))]
  vHalf = c(0, 0, 1, -1, rep(0, n-4))/sqrt(2)
  v0 = c(1, -1, rep(0, n-2))/sqrt(2)
  # -2 times the REML log likelihood, profiled and less a constant, from
  # the dense covariance matrix of the records, K/lambda + I
  X = matrix(1, n, 1)
  m2ll = function(lambda, y){
    Hi = solve(K/lambda + diag(n))
    XHX = t(X) %*% Hi %*% X
    Py = Hi %*% y - Hi %*% X %*% solve(XHX, t(X) %*% Hi %*% y)
    (n-1)*log(sum(y*Py)) -
      c(determinant(Hi)$modulus) + c(determinant(XHX)$modulus)
  }
  check = function(y, guess){
    ans = suppressWarnings(AlphaSimR:::.fitPedBLUP(y, iid, rep(1L, n), 2L,
                                                   simParam=SP,
                                                   lambdaGuess=guess))
    lambda = ans$Ve/ans$Vu
    grid = guess*10^seq(-2, 2, length.out=9)
    expect_gte(lambda, guess/100*(1-1e-6))
    expect_lte(lambda, guess*100*(1+1e-6))
    expect_lte(m2ll(lambda, y), min(vapply(grid, m2ll, numeric(1), y=y)) +
                 1e-6)
  }
  for(guess in c(1e-3, 1, 1e3)){
    check(sqrt(1e-6)*v0 + vHalf + sqrt(1e-6)*v4, guess)
    check(sqrt(3e-4)*v0 + vHalf + sqrt(1e-6)*v4, guess)
  }
})

test_that("GUESS the REML search is centered on the simulation's variance ratio", {
  sim = blupSim(9211)
  SP = sim$SP
  guess = AlphaSimR:::.lambdaGuess(1, "pheno", SP)
  expect_equal(guess, SP$varE[1]/SP$varA[1])
  expect_equal(guess, 0.6/0.4)
  expect_equal(AlphaSimR:::.lambdaGuess(SP$traitNames[1], "pheno", SP),
               guess)
  # Records that the trait's variances do not describe start from 1
  expect_equal(AlphaSimR:::.lambdaGuess(1, "gv", SP), 1)
  expect_equal(AlphaSimR:::.lambdaGuess(function(x) x[,1], "pheno", SP), 1)
  sol = pedBLUP(sim$train, simParam=SP)
  lambda = c(sol@Ve/sol@Vu)
  expect_true(lambda>=guess/100 && lambda<=guess*100)
})

test_that("LEGACY a SimParam that does not record resets is refused", {
  # A SimParam saved by an older version has no pedResets method, which
  # asking for it reports as NULL, as it does for this stand-in
  expect_error(AlphaSimR:::.pedResets(list()), "older version")
})

test_that("BOUNDARY a variance ratio at the edge of the range warns", {
  # Pairs of full sibs whose records are exact opposites have a negative
  # correlation, so the REML estimate of Vu is at zero
  set.seed(9106)
  founderPop = quickHaplo(nInd=20, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackPed(TRUE)
  SP$setSexes("no")
  pop = newPop(founderPop, simParam=SP)
  kids = makeCross(pop, cbind(seq(1, 19, 2), seq(2, 20, 2)), nProgeny=2,
                   simParam=SP)
  y = rep(c(1, -1), 10)
  expect_warning(ans <- AlphaSimR:::.fitPedBLUP(y, kids@iid,
                                                rep(1L, kids@nInd), 2L,
                                                simParam=SP),
                 "edge of the range")
  expect_true(ans$boundary)
})

# 5. pedBLUP and setEBV ----

test_that("EXPORT pedBLUP returns the internal fit as a PedSol", {
  sim = blupSim(9201)
  tr = sim$train
  sol = pedBLUP(tr, simParam=sim$SP)
  # pedBLUP starts the REML search from the simulation's variance ratio
  ref = fit(sim, lambdaGuess=AlphaSimR:::.lambdaGuess(1, "pheno", sim$SP))
  expect_s4_class(sol, "PedSol")
  expect_equal(c(sol@Vu), ref$Vu)
  expect_equal(c(sol@Ve), ref$Ve)
  expect_equal(sol@ebv, ref$ebv)
  expect_identical(sol@iid, as.integer(ref$iid))
  expect_equal(sol@center, mean(ref$ebv[match(tr@iid, ref$iid)]))
  expect_identical(sol@trait, sim$SP$traitNames[1])
})

test_that("SETEBV gv adds the first fixed effect and bv centers", {
  sim = blupSim(9202)
  tr = sim$train
  sol = pedBLUP(tr, simParam=sim$SP)
  u = sol@ebv[match(tr@iid, sol@iid)]
  gv = setEBV(tr, sol, value="gv", onlyValues=TRUE, simParam=sim$SP)
  bv = setEBV(tr, sol, value="bv", onlyValues=TRUE, simParam=sim$SP)
  expect_equal(c(gv), u + sol@fixEff[[1]])
  expect_equal(c(bv), u - sol@center)
  expect_equal(mean(bv), 0)
  expect_identical(colnames(gv), paste0("est_GV_", sol@trait))
  expect_identical(colnames(bv), paste0("est_BV_", sol@trait))
  # targetPop centers bv on its own mean
  last = tr[tr@fixEff==4L]
  bvT = setEBV(last, sol, value="bv", targetPop=last, onlyValues=TRUE,
               simParam=sim$SP)
  expect_equal(mean(bvT), 0)
  # The slot is replaced or appended as for an RRsol
  pop = setEBV(tr, sol, simParam=sim$SP)
  pop = setEBV(pop, sol, value="bv", append=TRUE, simParam=sim$SP)
  expect_identical(ncol(pop@ebv), 2L)
})

test_that("PREDICT new individuals get their BLUP from parent averages", {
  # Individuals created after the fit, including doubled haploids, crosses
  # of crosses and their progeny, against Henderson's BLUP with A extended
  # to them and no records of their own
  sim = blupSim(9203)
  tr = sim$train
  sol = pedBLUP(tr, Vu=0.7, Ve=1.3, simParam=sim$SP)
  last = tr[tr@fixEff==4L]
  new1 = makeDH(last[1:5], simParam=sim$SP)
  new2 = randCross(last, nCrosses=10, simParam=sim$SP)
  new3 = randCross(c(new1, new2), nCrosses=10, simParam=sim$SP)
  newPop = c(new1, new2, new3)
  pred = setEBV(newPop, sol, value="bv", onlyValues=TRUE, simParam=sim$SP)
  all = c(tr, newPop)
  ped = AlphaSimR:::.prunePed(all@iid, 2L, sim$SP)
  A = AlphaSimR:::calcPedA(ped$mother, ped$father, ped$isDH)
  X = fixedDesign(tr@fixEff)
  Z = matrix(0, tr@nInd, length(ped$iid))
  Z[cbind(seq_len(tr@nInd), ped$row[seq_len(tr@nInd)])] = 1
  ref = denseBLUP(tr@pheno[,1], X, Z, A, 0.7, 1.3)
  uNew = ref$u[ped$row[tr@nInd + seq_len(newPop@nInd)]]
  expect_equal(c(pred), uNew - sol@center, tolerance=1e-8)
})

test_that("VARCOMP an earlier fit's variance components are reused", {
  sim = blupSim(9204)
  tr = sim$train
  sol = pedBLUP(tr, simParam=sim$SP)
  sol2 = pedBLUP(tr, varComp=sol, simParam=sim$SP)
  expect_equal(c(sol2@Vu, sol2@Ve), c(sol@Vu, sol@Ve))
  expect_equal(sol2@ebv, sol@ebv, tolerance=1e-8)
})

test_that("SHIFT a large offset in the records changes only the fixed effects", {
  # Uncentered records sharing an offset of 2^30 lost their variation to
  # cancellation in y'Py
  sim = blupSim(9205)
  tr = sim$train
  y = tr@pheno[,1]
  f = function(y) AlphaSimR:::.fitPedBLUP(y, tr@iid, tr@fixEff, 2L,
                                          simParam=sim$SP)
  a = f(y)
  b = f(y + 2^30)
  expect_equal(b$Vu, a$Vu, tolerance=1e-5)
  expect_equal(b$Ve, a$Ve, tolerance=1e-5)
  expect_equal(b$ebv, a$ebv, tolerance=1e-5)
  expect_equal(unname(b$fixEff - 2^30), unname(a$fixEff), tolerance=1e-5)
  # Rescaling the records rescales the estimates and nothing else
  s = f(y*1e-8)
  expect_equal(s$Vu, a$Vu*1e-16, tolerance=1e-5)
  expect_equal(s$Ve, a$Ve*1e-16, tolerance=1e-5)
  expect_equal(s$ebv, a$ebv*1e-8, tolerance=1e-5)
})

test_that("RESET a fit is refused once a reset discards its individuals", {
  sim = blupSim(9206)
  tr = sim$train
  sol = pedBLUP(tr, Vu=1, Ve=1, simParam=sim$SP)
  # A reset that keeps every fitted individual leaves the fit usable
  sim$SP$resetPed(sim$SP$lastId)
  expect_silent(setEBV(tr, sol, simParam=sim$SP))
  # A reset below them lets new individuals take their iids
  founders = tr[tr@fixEff==1L][1:40]
  sim$SP$resetPed(max(founders@iid))
  newPop = randCross(founders, nCrosses=10, simParam=sim$SP)
  expect_true(all(newPop@iid %in% sol@iid))
  expect_error(setEBV(newPop, sol, simParam=sim$SP), "reset")
})

test_that("GIVEN variance components need no residual degrees of freedom", {
  # One unrelated individual with one record and an intercept: the
  # equations [1 1; 1 2] are nonsingular, giving u = 0 and b = y
  set.seed(9207)
  founderPop = quickHaplo(nInd=1, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackPed(TRUE)
  pop = newPop(founderPop, simParam=SP)
  ans = AlphaSimR:::.fitPedBLUP(5, pop@iid, 1L, 2L, Vu=1, Ve=1,
                                simParam=SP)
  expect_equal(ans$ebv, 0)
  expect_equal(unname(ans$fixEff), 5)
  expect_error(AlphaSimR:::.fitPedBLUP(5, pop@iid, 1L, 2L, simParam=SP),
               "more records")
  expect_error(AlphaSimR:::.fitPedBLUP(numeric(0), integer(0), integer(0),
                                       2L, simParam=SP),
               "no records")
})

test_that("CONSTANT records need given variance components", {
  # Records that do not vary within the levels of the fixed effects leave
  # nothing for REML to estimate, but fit with given components
  set.seed(9208)
  founderPop = quickHaplo(nInd=4, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackPed(TRUE)
  pop = newPop(founderPop, simParam=SP)
  two = pop@iid[1:2]
  expect_error(AlphaSimR:::.fitPedBLUP(c(5, 5), two, c(1L, 1L), 2L,
                                       simParam=SP),
               "do not vary")
  ans = AlphaSimR:::.fitPedBLUP(c(5, 5), two, c(1L, 1L), 2L, Vu=1,
                                Ve=1, simParam=SP)
  expect_equal(ans$ebv, c(0, 0))
  expect_equal(unname(ans$fixEff), 5)
  # Records that differ by far less than their size, or that are all far
  # below 1, still vary. With unrelated individuals the variance ratio is
  # not identified, so a warning about the range searched is expected.
  f = function(y) suppressWarnings(
    AlphaSimR:::.fitPedBLUP(y, pop@iid, rep(1L, 4), 2L, simParam=SP))
  expect_error(f(1e16 + c(-4, -2, 2, 4)), NA)
  expect_error(f(1e-16*c(-4, -2, 2, 4)), NA)
})

test_that("INBRED highly inbred selfed siblings keep their own breeding values", {
  # After 33 generations of selfing, a selfed offspring's Mendelian
  # sampling variance is 2^-34, which is genuine and must not be taken for
  # rounding error
  set.seed(9209)
  founderPop = quickHaplo(nInd=1, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackPed(TRUE)
  pop = newPop(founderPop, simParam=SP)
  for(g in 1:33){
    pop = self(pop, simParam=SP)
  }
  sibs = self(pop, nProgeny=2, simParam=SP)
  ped = SP$pedigree
  ans = AlphaSimR:::calcPedAinvTriplet(ped[,1], ped[,2], ped[,3])
  expect_true(all(sibs@iid %in% ans$basis))
  expect_equal(ans$D[sibs@iid], rep(2^-34, 2))
  # With the variance treated as zero, the siblings' breeding values were
  # forced to be equal, whatever Ve. With Ve comparable to it, each takes
  # about D/(D+Ve) of its deviation from their mean.
  y = c(1, -1)
  fitted = AlphaSimR:::.fitPedBLUP(y, sibs@iid, c(1L, 1L), 2L, Vu=1,
                                   Ve=1e-10, simParam=SP)
  pp = AlphaSimR:::.prunePed(sibs@iid, 2L, SP)
  A = AlphaSimR:::calcPedA(pp$mother, pp$father, pp$isDH)
  Z = matrix(0, 2, length(pp$iid))
  Z[cbind(1:2, pp$row)] = 1
  ref = denseBLUP(y, matrix(1, 2, 1), Z, A, 1, 1e-10)
  u = fitted$ebv[pp$row]
  # The equations are badly conditioned here, and both solutions differ
  # by about 1e-6 between platforms
  expect_equal(u, ref$u[pp$row], tolerance=1e-5)
  expect_equal(u[1], 2^-34/(2^-34+1e-10), tolerance=1e-4)
  expect_equal(u[1], -u[2], tolerance=1e-5)
  # Twelve more generations bring the variance to 2^-46, still genuine.
  # With Ve equal to it, each sibling takes half of its deviation from
  # their mean. A dense reference cannot resolve differences this small.
  for(g in 1:12){
    pop = self(pop, simParam=SP)
  }
  sibs = self(pop, nProgeny=2, simParam=SP)
  ped = SP$pedigree
  ans = AlphaSimR:::calcPedAinvTriplet(ped[,1], ped[,2], ped[,3])
  expect_true(all(sibs@iid %in% ans$basis))
  expect_equal(ans$D[sibs@iid], rep(2^-46, 2))
  fitted = AlphaSimR:::.fitPedBLUP(y, sibs@iid, c(1L, 1L), 2L, Vu=1,
                                   Ve=2^-46, simParam=SP)
  pp = AlphaSimR:::.prunePed(sibs@iid, 2L, SP)
  expect_equal(fitted$ebv[pp$row], c(0.5, -0.5), tolerance=1e-4)
  # At 2^-53 the variances along the selfing chain span the whole precision
  # of a double, and whether the equations factor depends on the platform's
  # rounding. REML leaves out the ratios at which they cannot be factored,
  # such as a ratio of 1, where the records round away, and if none can be,
  # it says so rather than stopping with CHOLMOD's message.
  for(g in 1:7){
    pop = self(pop, simParam=SP)
  }
  sibs = self(pop, nProgeny=2, simParam=SP)
  ped = SP$pedigree
  ans = AlphaSimR:::calcPedAinvTriplet(ped[,1], ped[,2], ped[,3])
  expect_identical(c(ans$D)[sibs@iid], rep(2^-53, 2))
  expect_true(all(sibs@iid %in% ans$basis))
  res = tryCatch(suppressWarnings(
    AlphaSimR:::.fitPedBLUP(y, sibs@iid, c(1L, 1L), 2L, simParam=SP)),
    error=function(e) e)
  if(inherits(res, "error")){
    expect_match(conditionMessage(res), "too badly conditioned")
  }else{
    expect_true(is.finite(res$Vu) && is.finite(res$Ve))
    expect_true(all(is.finite(res$ebv)))
  }
  # A cross of a parent selfed for 53 generations, F = 1-2^-53, with a
  # fully inbred one has a variance of 2^-55, which the sum of their
  # inbreeding coefficients, 2-2^-53, once rounded to zero
  pop = self(pop, simParam=SP)
  dh = makeDH(pop, simParam=SP)
  sibs = makeCross(c(pop, dh), crossPlan=cbind(1L, 2L), nProgeny=2,
                   simParam=SP)
  ped = SP$pedigree
  ans = AlphaSimR:::calcPedAinvTriplet(ped[,1], ped[,2], ped[,3])
  # Exact comparisons, as the default tolerance cannot tell these from 1
  # and 0
  expect_identical(c(ans$F)[pop@iid], 1-2^-53)
  expect_identical(c(ans$D)[sibs@iid], rep(2^-55, 2))
  expect_true(all(sibs@iid %in% ans$basis))
  # BLUP is not checked here, as equations this badly conditioned factor on
  # some platforms and not on others; the variance and the basis are exact
})

test_that("INBRED clones of inbred lines still have no Mendelian sampling", {
  # Selfs and doubled haploids of an inbred founder are copies of it, and
  # an F1 of two inbred founders is their average, whatever the depth
  mother = c(0L, 0L, 1L, 3L, 1L, 4L, 6L)
  father = c(0L, 0L, 1L, 3L, 2L, 4L, 5L)
  isDH =   c(1L, 1L, 0L, 1L, 0L, 0L, 0L)
  ans = AlphaSimR:::calcPedAinvTriplet(mother, father, isDH)
  # 3 is a self of 1, 4 a doubled haploid of 3 and 6 a self of 4, all
  # copies of 1; 5 is the F1 of 1 and 2; 7 crosses a copy of 1 with 5
  expect_equal(ans$basis, c(1L, 2L, 7L))
  expect_equal(c(ans$D)[3:6], c(0, 0, 0, 0))
  expect_equal(c(ans$F)[3:6], c(1, 1, 0, 1))
})

# 6. Errors ----

test_that("ERROR inputs are checked", {
  sim = blupSim(9107)
  tr = sim$train
  y = tr@pheno[,1]
  f = function(...) AlphaSimR:::.fitPedBLUP(..., simParam=sim$SP)
  yNA = y
  yNA[3] = NA
  expect_error(f(yNA, tr@iid, tr@fixEff, 2L), "missing values")
  expect_error(f(y, tr@iid, tr@fixEff, 2L, Vu=1), "both Vu and Ve")
  expect_error(f(y, tr@iid, tr@fixEff, 2L, Vu=-1, Ve=1),
               "single positive value")
  expect_error(f(y[1:3], tr@iid[1:3], 1:3, 2L), "more records")
  expect_error(f(y, tr@iid[-1], tr@fixEff, 2L), "same length")
  sol = pedBLUP(tr, simParam=sim$SP)
  expect_error(pedBLUP(tr, traits=1:2, simParam=sim$SP), "one trait")
  expect_error(pedBLUP(tr, Vu=1, varComp=sol, simParam=sim$SP), "not both")
  expect_error(pedBLUP(tr, varComp=list(), simParam=sim$SP), "earlier fit")
  expect_error(setEBV(tr, sol, value="female", simParam=sim$SP), "GCA")
})
