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
  expect_equal(ans$Vu, c(ref$Vu), tolerance=1e-4)
  expect_equal(ans$Ve, c(ref$Ve), tolerance=1e-4)
  expect_false(ans$boundary)
  # The BLUP at the estimates matches the dense BLUP at them
  d = denseDesign(sim)
  den = denseBLUP(tr@pheno[,1], d$X, d$Z, d$A, ans$Vu, ans$Ve)
  expect_equal(ans$ebv, den$u, tolerance=1e-8)
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

# 4. Errors ----

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
})
