context("GS")

# The genomic selection functions had no tests at all. These do not check the
# models against an outside implementation, because there isn't one here.
# They use the redundancy inside the module instead: the same model can be
# reached by more than one route, and the routes have to agree.
#
# The most important of those is RRBLUP given its variance components. It
# solves the mixed model in one of two ways, working with a square matrix of
# the records when they are fewer than the markers and of the markers
# otherwise, and the pair only means anything if both give the same answer as
# the REML fit does.

# The numbered models are deprecated, and are still tested until they are
# removed, without the warning they give
quietDeprecated = function(expr){
  suppressWarnings(expr, classes="deprecatedWarning")
}

# A population with more records than markers, so that the fixed variance
# solve works with the markers, and one with fewer, so that it works with the
# records
gsPop = function(nInd, nQtl, nSnp, nChr=2, segSites=NULL, seed=8001,
                 nFixEff=1L){
  if(is.null(segSites)){
    segSites = nQtl + nSnp
  }
  set.seed(seed)
  founderPop = quickHaplo(nInd=nInd, nChr=nChr, segSites=segSites)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$restrSegSites(minQtlPerChr=nQtl, minSnpPerChr=nSnp, overlap=FALSE)
  SP$addTraitA(nQtlPerChr=nQtl)
  SP$addSnpChip(nSnpPerChr=nSnp)
  SP$setVarE(h2=0.5)
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)
  if(nFixEff > 1L){
    pop@fixEff = rep_len(seq_len(nFixEff), pop@nInd)
  }
  return(list(pop=pop, SP=SP, nMarker=nChr*nSnp))
}

# Marker effects from a solution, as a plain vector
addEff = function(ans, i=1L){
  return(c(ans@gv[[i]]@addEff))
}

test_that("RRBLUP given its variance components agrees with REML when it works with the markers", {
  # 40 records against 10 markers, so Henderson's equations are cheaper
  d = gsPop(nInd=40, nQtl=5, nSnp=5)
  expect_lt(d$nMarker, nInd(d$pop))

  ref = RRBLUP(d$pop, simParam=d$SP)
  ans = RRBLUP(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1],
               simParam=d$SP)

  expect_equal(addEff(ans), addEff(ref), tolerance=1e-6)
  expect_equal(ans@gv[[1]]@intercept, ref@gv[[1]]@intercept, tolerance=1e-6)
})

test_that("RRBLUP given its variance components agrees with REML when it works with the records", {
  # 15 records against 40 markers, so the record side is cheaper. This is the
  # branch that never forms a markers by records matrix, and the one a
  # mistake in the order of that product would show up in.
  d = gsPop(nInd=15, nQtl=5, nSnp=20)
  expect_gt(d$nMarker, nInd(d$pop))

  ref = RRBLUP(d$pop, simParam=d$SP)
  ans = RRBLUP(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1],
               simParam=d$SP)

  expect_equal(addEff(ans), addEff(ref), tolerance=1e-6)
  expect_equal(ans@gv[[1]]@intercept, ref@gv[[1]]@intercept, tolerance=1e-6)
})

test_that("both sides of the fixed variance solve handle more than one fixed effect level", {
  # The record side splits the solution of V into the fixed effect columns
  # and the response, so it has to keep them apart when there is more than
  # one fixed effect. Both shapes of the problem are checked.
  for(size in list(c(nInd=40, nQtl=5, nSnp=5),
                   c(nInd=15, nQtl=5, nSnp=20))){
    d = gsPop(nInd=size[["nInd"]], nQtl=size[["nQtl"]], nSnp=size[["nSnp"]],
              nFixEff=3L)
    expect_equal(length(unique(d$pop@fixEff)), 3L)

    ref = RRBLUP(d$pop, simParam=d$SP)
    ans = RRBLUP(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1],
                 simParam=d$SP)

    expect_equal(addEff(ans), addEff(ref), tolerance=1e-6,
                 info=paste("nInd", size[["nInd"]]))
    expect_true(all(is.finite(addEff(ans))))
  }
})

test_that("RRBLUP2 estimates variance components when asked to", {
  # The EM path is the only one that asks the factorization for a trace, so
  # it exercises code the other path never reaches. There is nothing to
  # compare the estimates against, so this checks that they are usable.
  d = gsPop(nInd=40, nQtl=5, nSnp=5)

  ans = quietDeprecated(RRBLUP2(d$pop, useEM=TRUE, maxIter=200,
                                simParam=d$SP))

  expect_true(is.finite(ans@Vu[1,1]))
  expect_true(is.finite(ans@Ve[1,1]))
  expect_gt(ans@Vu[1,1], 0)
  expect_gt(ans@Ve[1,1], 0)
  expect_true(all(is.finite(addEff(ans))))

  # Predictions from an estimated model should still track the truth
  pred = setEBV(d$pop, ans, simParam=d$SP)
  expect_gt(cor(c(ebv(pred)), c(gv(d$pop))), 0.3)
})

test_that("RRBLUP given its variance components holds them fixed", {
  # Supplying the variance components turns estimation off, so two calls with
  # the same components must agree exactly and a third with different
  # components must not
  d = gsPop(nInd=30, nQtl=5, nSnp=5)
  ref = RRBLUP(d$pop, simParam=d$SP)

  a = RRBLUP(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1],
             simParam=d$SP)
  b = RRBLUP(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1],
             simParam=d$SP)
  expect_identical(addEff(a), addEff(b))

  c2 = RRBLUP(d$pop, Vu=ref@Vu[1,1]*10, Ve=ref@Ve[1,1],
              simParam=d$SP)
  expect_false(isTRUE(all.equal(addEff(a), addEff(c2))))

  # A larger marker variance shrinks the effects less
  expect_gt(sum(addEff(c2)^2), sum(addEff(a)^2))
})

test_that("fast mode RRBLUP agrees with RRBLUP given the same variance components", {
  # Iterates to convergence rather than solving, so it is the slowest fit here
  skip_on_cran()
  d = gsPop(nInd=40, nQtl=5, nSnp=5)
  ref = RRBLUP(d$pop, simParam=d$SP)
  ans = RRBLUP(d$pop, fastMode=TRUE, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1],
               maxIter=10000, simParam=d$SP)

  # Fast mode iterates rather than solving directly, so it is held to a
  # looser tolerance than the direct solvers are held to each other
  expect_equal(addEff(ans), addEff(ref), tolerance=1e-4)
})

# Single trait fast mode deflates the largest eigenvalues of its Jacobi
# preconditioned system, using eigenvectors estimated from a subset of the
# records. The preconditioner changes how fast the iteration converges, never
# what it converges to, so the solution has to be the same with deflation as
# without it. The exported solver is called directly here, because it
# reports how many eigenvalues it deflated and how many iterations it took.
fastSolve = function(d, Vu, Ve, nPrecond, maxRank, nThreads=1L){
  snp = d$SP$snpChips[[1]]
  return(AlphaSimR:::callFastRRBLUP(d$pop@pheno, rep(1L, nInd(d$pop)), d$pop,
                                    snp@lociPerChr, snp@lociLoc, Vu, Ve,
                                    10000L, FALSE, 1L, nPrecond, maxRank,
                                    nThreads))
}

# 200 records against 100 markers in ten full-sib families. Family structure
# is what puts eigenvalues of the scaled equations well above twice their
# mean, so that some are deflated; unrelated founders from quickHaplo have
# hardly any.
famPop = function(){
  d = gsPop(nInd=40, nQtl=10, nSnp=50)
  d$pop = randCross(d$pop, nCrosses=10, nProgeny=20, simParam=d$SP)
  d$pop = setPheno(d$pop, simParam=d$SP)
  return(d)
}

test_that("fast mode RRBLUP converges to the same solution with deflation", {
  skip_on_cran()
  d = famPop()
  Vu = 0.1
  Ve = 1
  ref = RRBLUP(d$pop, Vu=Vu, Ve=Ve, simParam=d$SP)

  jacobi = fastSolve(d, Vu, Ve, 0L, 0L)
  expect_equal(jacobi$rank, 0L)
  # The whole population and a subset of it, which estimates the
  # eigenvectors less well
  for(nPrecond in c(200L, 60L)){
    deflated = fastSolve(d, Vu, Ve, nPrecond, 1000L)
    expect_gt(deflated$rank, 0L)
    expect_lt(deflated$rank, nPrecond)
    expect_lte(deflated$iter, jacobi$iter)
    expect_equal(c(deflated$alpha), c(jacobi$alpha), tolerance=1e-6,
                 info=paste("nPrecond", nPrecond))
    expect_equal(c(deflated$alpha), addEff(ref), tolerance=1e-4,
                 info=paste("nPrecond", nPrecond))
  }

  # The cap on the rank is respected
  capped = fastSolve(d, Vu, Ve, 200L, 3L)
  expect_equal(capped$rank, 3L)
  expect_equal(c(capped$alpha), c(jacobi$alpha), tolerance=1e-6)
})

test_that("fast mode RRBLUP options turn deflation off and cap its rank", {
  # Both options are saved, so that neither is left set for later tests
  old = options(AlphaSimR.fastModeDeflation=FALSE,
                AlphaSimR.fastModeMaxRank=NULL)
  on.exit(options(old))
  expect_equal(AlphaSimR:::fastDeflation(3000L, 100L)$maxRank, 0L)
  options(AlphaSimR.fastModeDeflation=NULL, AlphaSimR.fastModeMaxRank=7L)
  expect_equal(AlphaSimR:::fastDeflation(3000L, 100L)$maxRank, 7L)
  # Five times the square root of the number of records, but never more than
  # nSubsample, or 3000 when it is zero
  expect_equal(AlphaSimR:::fastDeflation(3000L, 100L)$nRecords, 50L)
  expect_equal(AlphaSimR:::fastDeflation(3000L, 1e5L)$nRecords, 1581L)
  expect_equal(AlphaSimR:::fastDeflation(0L, 1e6L)$nRecords, 3000L)
  expect_equal(AlphaSimR:::fastDeflation(500L, 1e5L)$nRecords, 500L)
})

test_that("fast mode RRBLUP given its variance components draws no random numbers", {
  # The records of the preconditioner are evenly spaced rather than sampled,
  # so a fit with supplied components leaves the random number stream where
  # it was
  d = famPop()
  set.seed(5)
  before = .Random.seed
  ans = RRBLUP(d$pop, fastMode=TRUE, Vu=0.01, Ve=1, simParam=d$SP)
  expect_identical(.Random.seed, before)
})

test_that("fast mode RRBLUP with deflation is not affected by the number of threads", {
  skip_on_cran()
  nThreads = getNumThreads()
  skip_if_not(nThreads > 1L, "only one thread available")
  d = famPop()
  one = fastSolve(d, 0.1, 1, 200L, 1000L, nThreads=1L)
  many = fastSolve(d, 0.1, 1, 200L, 1000L, nThreads=nThreads)
  expect_gt(one$rank, 0L)
  expect_identical(one$alpha, many$alpha)
  expect_identical(one$iter, many$iter)
})

test_that("RRBLUP outside fast mode takes both variance components or neither", {
  # Estimating the other one would ignore the one supplied without notice
  d = gsPop(nInd=30, nQtl=5, nSnp=5)
  expect_error(RRBLUP(d$pop, Vu=0.01, simParam=d$SP), "both Vu and Ve")
  expect_error(RRBLUP(d$pop, Ve=1, simParam=d$SP), "both Vu and Ve")
  expect_error(RRBLUP(d$pop, Vu=-0.01, Ve=1, simParam=d$SP),
               "must not be negative")
  expect_error(RRBLUP(d$pop, Vu=0.01, Ve=0, simParam=d$SP),
               "Ve must be positive")
})

test_that("setEBV reproduces the model's own predictions", {
  d = gsPop(nInd=30, nQtl=5, nSnp=5)
  ans = RRBLUP(d$pop, simParam=d$SP)
  pred = setEBV(d$pop, ans, simParam=d$SP)

  expect_equal(nrow(ebv(pred)), nInd(d$pop))
  expect_equal(ncol(ebv(pred)), 1L)
  expect_true(all(is.finite(c(ebv(pred)))))

  # Worked out by hand from the marker effects and the dosages. A locus
  # contributes its effect times (dosage - ploidy/2) * (2/ploidy), so the
  # dosages are centered before they are weighted rather than used raw.
  M = pullSnpGeno(d$pop, simParam=d$SP)
  a = addEff(ans)
  p = d$pop@ploidy
  centered = (M - p/2) * (2/p)
  byHand = c(centered %*% a) + ans@gv[[1]]@intercept
  expect_equal(unname(c(ebv(pred))), unname(byHand), tolerance=1e-6)

  # The model was fitted on these records, so it should fit them well
  expect_gt(cor(c(ebv(pred)), c(pheno(d$pop))), 0.5)
})

test_that("setEBV appends rather than replaces when asked", {
  d = gsPop(nInd=20, nQtl=5, nSnp=5)
  ans = RRBLUP(d$pop, simParam=d$SP)

  once = setEBV(d$pop, ans, simParam=d$SP)
  twice = setEBV(once, ans, append=TRUE, simParam=d$SP)

  expect_equal(ncol(ebv(once)), 1L)
  expect_equal(ncol(ebv(twice)), 2L)
  expect_equal(unname(ebv(twice)[,1]), unname(ebv(twice)[,2]))

  # Without append the second call replaces the first
  replaced = setEBV(once, ans, simParam=d$SP)
  expect_equal(ncol(ebv(replaced)), 1L)
})

test_that("setEBV returns only the estimates when onlyValues is TRUE", {
  d = gsPop(nInd=20, nQtl=5, nSnp=5)
  ans = RRBLUP(d$pop, simParam=d$SP)

  pred = setEBV(d$pop, ans, simParam=d$SP)
  est = setEBV(d$pop, ans, onlyValues=TRUE, simParam=d$SP)

  # The matrix is the one that would have been written to the EBV slot,
  # and append has no effect because nothing is written
  expect_true(is.matrix(est))
  expect_equal(est, ebv(pred))
  expect_equal(setEBV(pred, ans, append=TRUE, onlyValues=TRUE,
                      simParam=d$SP), est)
})

test_that("a model fitted on genetic values recovers the QTL effects", {
  # Needs a few hundred individuals for the effects to be recoverable
  skip_on_cran()
  # With no environmental noise and the QTL themselves as markers, the fitted
  # effects should line up with the true ones. This is the one check that
  # looks at whether the model is right rather than whether two routes agree.
  set.seed(8101)
  founderPop = quickHaplo(nInd=200, nChr=2, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=10)
  SP$setVarE(h2=0.999)
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)

  ans = RRBLUP(pop, useQtl=TRUE, simParam=SP)
  trueEff = SP$traits[[1]]@addEff

  expect_equal(length(addEff(ans)), length(trueEff))
  expect_gt(cor(addEff(ans), trueEff), 0.9)

  pred = setEBV(pop, ans, simParam=SP)
  expect_gt(cor(c(ebv(pred)), c(gv(pop))), 0.95)
})

test_that("the dominance and hybrid models run and predict", {
  # Three multi kernel fits, each holding a matrix of several random effects
  skip_on_cran()
  # These fit more than one random effect, which is a different solver. No
  # reference to compare against, so this checks they produce usable output.
  set.seed(8201)
  founderPop = quickHaplo(nInd=60, nChr=2, segSites=20)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setSexes("yes_sys")
  SP$restrSegSites(minQtlPerChr=10, minSnpPerChr=10, overlap=FALSE)
  SP$addTraitAD(nQtlPerChr=10, meanDD=0.5)
  SP$addSnpChip(nSnpPerChr=10)
  SP$setVarE(h2=0.5)
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)

  for(fit in list(RRBLUP_D, RRBLUP_GCA, RRBLUP_SCA)){
    ans = fit(pop, simParam=SP)
    expect_true(is(ans, "RRsol"))
    expect_true(all(is.finite(addEff(ans))))
    pred = setEBV(pop, ans, simParam=SP)
    expect_equal(nrow(ebv(pred)), nInd(pop))
    expect_true(all(is.finite(c(ebv(pred)))))
  }
})

test_that("the dominance and hybrid models run with fixed components", {
  # Three more multi kernel fits
  skip_on_cran()
  # useEM=FALSE is the path that no longer inverts its Cholesky factor
  set.seed(8301)
  founderPop = quickHaplo(nInd=60, nChr=2, segSites=20)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setSexes("yes_sys")
  SP$restrSegSites(minQtlPerChr=10, minSnpPerChr=10, overlap=FALSE)
  SP$addTraitAD(nQtlPerChr=10, meanDD=0.5)
  SP$addSnpChip(nSnpPerChr=10)
  SP$setVarE(h2=0.5)
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)

  vg = varG(pop)[1,1]
  ve = varP(pop)[1,1] - vg

  dom = RRBLUP_D(pop, Vu=c(vg/20, vg/20), Ve=ve, simParam=SP)
  expect_true(all(is.finite(addEff(dom))))

  gca = RRBLUP_GCA(pop, Vu=c(vg/20, vg/20), Ve=ve,
                   simParam=SP)
  expect_true(all(is.finite(addEff(gca))))

  sca = RRBLUP_SCA(pop, Vu=c(vg/20, vg/20, vg/20), Ve=ve,
                   simParam=SP)
  expect_true(all(is.finite(addEff(sca))))
})

test_that("a model can be fitted to more than one trait at once", {
  # Fits a multivariate model, which is the most expensive solver here
  skip_on_cran()
  set.seed(8401)
  founderPop = quickHaplo(nInd=50, nChr=2, segSites=20)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$restrSegSites(minQtlPerChr=10, minSnpPerChr=10, overlap=FALSE)
  SP$addTraitA(nQtlPerChr=10, mean=c(0,0), var=c(1,1),
               corA=matrix(c(1,0.5,0.5,1), nrow=2))
  SP$addSnpChip(nSnpPerChr=10)
  SP$setVarE(h2=c(0.5,0.5))
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)

  ans = RRBLUP(pop, traits=1:2, simParam=SP)
  expect_equal(length(ans@gv), 2L)

  pred = setEBV(pop, ans, simParam=SP)
  expect_equal(ncol(ebv(pred)), 2L)
  expect_true(all(is.finite(c(ebv(pred)))))

  # Each estimate should track its own trait more closely than the other
  expect_gt(cor(ebv(pred)[,1], gv(pop)[,1]),
            cor(ebv(pred)[,1], gv(pop)[,2]))
})

test_that("RRBLUPMemUse answers for every model it accepts", {
  models = c("fastRRBLUP","RRBLUP","RRBLUP2","RRBLUP_D","RRBLUP_D2",
             "RRBLUP_GCA","RRBLUP_GCA2","RRBLUP_SCA","RRBLUP_SCA2")
  for(m in models){
    v = RRBLUPMemUse(nInd=1000, nMarker=500, model=m)
    expect_true(is.finite(v), info=m)
    expect_gt(v, 0)
  }

  # The names used before the models were named after their functions
  expect_equal(RRBLUPMemUse(nInd=100, nMarker=50, model="REG"),
               RRBLUPMemUse(nInd=100, nMarker=50, model="RRBLUP"))
  expect_equal(RRBLUPMemUse(nInd=100, nMarker=50, model="GCA"),
               RRBLUPMemUse(nInd=100, nMarker=50, model="RRBLUP_GCA"))
  expect_equal(RRBLUPMemUse(nInd=100, nMarker=50, model="SCA"),
               RRBLUPMemUse(nInd=100, nMarker=50, model="RRBLUP_SCA"))
  # fastRRBLUP became RRBLUP's fast mode
  expect_equal(RRBLUPMemUse(nInd=100, nMarker=50, model="fastRRBLUP"),
               RRBLUPMemUse(nInd=100, nMarker=50, model="RRBLUP",
                            fastMode=TRUE))
  expect_false(isTRUE(all.equal(
    RRBLUPMemUse(nInd=1000, nMarker=500, model="RRBLUP"),
    RRBLUPMemUse(nInd=1000, nMarker=500, model="RRBLUP", fastMode=TRUE))))

  expect_error(RRBLUPMemUse(nInd=100, nMarker=50, model="notAModel"),
               "not recognized")

  # Every estimate grows with the size of the problem
  for(m in models){
    small = RRBLUPMemUse(nInd=500, nMarker=500, model=m)
    large = RRBLUPMemUse(nInd=2000, nMarker=2000, model=m)
    expect_gt(large, small)
  }
})

test_that("RRBLUPMemUse counts the deflation of single trait fast mode", {
  # The eigenvectors of the markers the deflation keeps are held through the
  # whole iteration, 1,000 of them by default for 100,000 records, which at
  # 50,000 markers is 0.4 GB
  on = RRBLUPMemUse(nInd=1e5, nMarker=5e4, fastMode=TRUE)
  old = options(AlphaSimR.fastModeDeflation=FALSE)
  on.exit(options(old))
  off = RRBLUPMemUse(nInd=1e5, nMarker=5e4, fastMode=TRUE)
  expect_equal(on-off, 5e4*1000*8/1e9, tolerance=1e-8)
})

test_that("RRBLUPMemUse counts the multiple trait optimizer", {
  # The Newton steps hold several square matrices over the t*(t+1)
  # parameters of Vu and Ve, so with many traits they outweigh a small data
  # set. Four of them at 50 traits are already 0.21 GB.
  nPar = 50*51
  floorGB = 4*8*nPar^2/10^9
  for(m in c(50, 500)){
    expect_gt(RRBLUPMemUse(nInd=200, nMarker=m, nTraits=50), floorGB)
  }
})

test_that("RRBLUPMemUse counts the multiple trait fixed effects", {
  # Each likelihood evaluation holds a q by q inverse per trait, and the
  # current and trial points coexist, so 500 fixed effects on two traits
  # need at least 8 MB on either decomposition
  floorGB = 2*2*8*500^2/10^9
  for(m in c(10, 2000)){
    few = RRBLUPMemUse(nInd=1000, nMarker=m, nTraits=2, nFixEff=1)
    many = RRBLUPMemUse(nInd=1000, nMarker=m, nTraits=2, nFixEff=500)
    expect_gt(many - few, floorGB)
  }
})

test_that("RRBLUPMemUse counts the multiple trait multikernel models", {
  # With 1000 records and 200 markers per term, three traits put the SCA
  # model on Henderson's equations with 3*(1+3*200) coefficients, whose
  # coefficient matrix and its factors are held five times over
  nCoef = 3*(1+3*200)
  floorGB = 5*8*nCoef^2/10^9
  three = RRBLUPMemUse(nInd=1000, nMarker=200, model="RRBLUP_SCA", nTraits=3)
  one = RRBLUPMemUse(nInd=1000, nMarker=200, model="RRBLUP_SCA", nTraits=1)
  expect_gt(three, floorGB)
  expect_gt(three, one)
  for(m in c("RRBLUP_D","RRBLUP_GCA","fastRRBLUP")){
    expect_gt(RRBLUPMemUse(nInd=1000, nMarker=200, model=m, nTraits=3),
              RRBLUPMemUse(nInd=1000, nMarker=200, model=m, nTraits=1),
              label=m)
  }
})

test_that("RRBLUPMemUse counts the multiple trait working variates", {
  # The derivatives hold a vector of records per parameter twice over. For
  # the SCA model with 10 traits on 10,000 records that is two matrices of
  # 100,000 rows by 4*55 parameters.
  floorGB = 2*8*(10000*10)*(4*55)/10^9
  expect_gt(RRBLUPMemUse(nInd=10000, nMarker=10, model="RRBLUP_SCA",
                         nTraits=10), floorGB)
})

test_that("RRBLUPMemUse follows the multiple trait models to the record side", {
  # With many more markers than records, the record side is the cheaper
  # one, and its square matrix has a row per record and trait
  nRec = 2*200
  floorGB = 5*8*nRec^2/10^9
  v = RRBLUPMemUse(nInd=200, nMarker=5000, model="RRBLUP_D", nTraits=2)
  expect_gt(v, floorGB)
  expect_lt(v, RRBLUPMemUse(nInd=200, nMarker=5000, model="RRBLUP_D2",
                            nTraits=1))
  # And fast mode RRBLUP's variance components from a record cross product
  expect_gt(RRBLUPMemUse(nInd=10000, nMarker=10000, model="RRBLUP",
                         fastMode=TRUE, nTraits=3, nSubsample=2000),
            4*8*2000^2/10^9)
})

test_that("RRBLUPMemUse follows RRBLUP2 between its two methods", {
  # Estimating variance components holds one more square matrix
  withEM = RRBLUPMemUse(nInd=1000, nMarker=500, model="RRBLUP2", useEM=TRUE)
  noEM = RRBLUPMemUse(nInd=1000, nMarker=500, model="RRBLUP2", useEM=FALSE)
  expect_gt(withEM, noEM)

  # With more markers than records and the components fixed, the estimate is
  # for the record side and so does not grow with the markers the way the
  # marker side does
  recordSide = RRBLUPMemUse(nInd=500, nMarker=5000, model="RRBLUP2",
                            useEM=FALSE)
  markerSide = RRBLUPMemUse(nInd=500, nMarker=5000, model="RRBLUP2",
                            useEM=TRUE)
  expect_lt(recordSide, markerSide)

  # RRBLUP decomposes whichever cross product is smaller, so its estimate is
  # near the smaller of the two dimensions
  expect_lt(RRBLUPMemUse(nInd=500, nMarker=5000, model="RRBLUP"),
            RRBLUPMemUse(nInd=5000, nMarker=5000, model="RRBLUP"))
})

test_that("RRBLUPMemUse follows the multikernel models to Henderson's equations", {
  # With far more records than marker effects the dominance, GCA and SCA
  # models hold no square matrix of the records, so doubling the records
  # roughly doubles the estimate rather than quadrupling it
  for(m in c("RRBLUP_D","RRBLUP_GCA","RRBLUP_SCA")){
    small = RRBLUPMemUse(nInd=20000, nMarker=500, model=m)
    large = RRBLUPMemUse(nInd=40000, nMarker=500, model=m)
    expect_lt(large/small, 2.5, label=m)
  }
  # With more marker effects than records they hold square matrices of the
  # records, which grow with the square of the records
  small = RRBLUPMemUse(nInd=1000, nMarker=5000, model="RRBLUP_SCA")
  large = RRBLUPMemUse(nInd=2000, nMarker=5000, model="RRBLUP_SCA")
  expect_gt(large/small, 2.5)
})

test_that("models can be fitted on values other than phenotypes", {
  d = gsPop(nInd=40, nQtl=5, nSnp=5)
  for(use in c("pheno","gv","rand")){
    ans = RRBLUP(d$pop, use=use, simParam=d$SP)
    expect_true(all(is.finite(addEff(ans))), info=use)
  }

  # Whatever it is asked for comes back as a one column matrix, which is
  # what the models read the number of traits from
  for(use in c("pheno","gv","rand")){
    y = getResponse(pop=d$pop, trait=1, use=use, simParam=d$SP)
    expect_true(is.matrix(y), info=use)
    expect_equal(ncol(y), 1L, info=use)
    expect_equal(nrow(y), nInd(d$pop), info=use)
  }

  # A model fitted on genetic values predicts them better than one fitted on
  # phenotypes, because there is no noise in the response
  onGv = setEBV(d$pop, RRBLUP(d$pop, use="gv", simParam=d$SP), simParam=d$SP)
  onPheno = setEBV(d$pop, RRBLUP(d$pop, use="pheno", simParam=d$SP),
                   simParam=d$SP)
  expect_gt(cor(c(ebv(onGv)), c(gv(d$pop))),
            cor(c(ebv(onPheno)), c(gv(d$pop))))
})

test_that("fitting is not affected by the number of threads", {
  skip_on_cran()
  nThreads = getNumThreads()
  skip_if_not(nThreads > 1L, "only one thread available")

  d = gsPop(nInd=40, nQtl=5, nSnp=5)
  one = RRBLUP(d$pop, nThreads=1L, simParam=d$SP)
  many = RRBLUP(d$pop, nThreads=nThreads, simParam=d$SP)
  expect_equal(addEff(one), addEff(many), tolerance=1e-8)

  ref = one
  oneB = RRBLUP(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1],
                nThreads=1L, simParam=d$SP)
  manyB = RRBLUP(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1],
                 nThreads=nThreads, simParam=d$SP)
  expect_equal(addEff(oneB), addEff(manyB), tolerance=1e-8)
})

# The models read their markers straight from the population's packed
# genotypes. Building the same marker matrices in R from the pull functions
# and handing them to the exported solvers is an independent route to the
# same fit, and with a single fixed effect level the design matrix is a
# column of ones on both routes. The codings below are exact in double
# precision, so a fit through the same solver agrees to the last bit. A fit
# given its variance components is held to rounding against the exported
# solver for them, which is a different one.

# Dominance coding of a dosage x at ploidy P, as the models code it
domCode = function(x, P) x*(P - x)*(2/P)^2

test_that("RRBLUP reads the same markers that pullSnpGeno returns", {
  d = gsPop(nInd=60, nQtl=8, nSnp=12, nChr=3)
  y = d$pop@pheno[, 1L, drop=FALSE]
  M = pullSnpGeno(d$pop, simParam=d$SP) - 1
  ref = solveRRBLUP(y, matrix(1, nrow(y), 1L), M)
  ans = RRBLUP(d$pop, simParam=d$SP)
  expect_identical(addEff(ans), c(ref$u))

  ref2 = solveRRBLUP_EM(y, matrix(1, nrow(y), 1L), M, Vu=0.1, Ve=1,
                        tol=1e-6, maxIter=10L, useEM=FALSE)
  ans2 = RRBLUP(d$pop, Vu=0.1, Ve=1, simParam=d$SP)
  expect_equal(addEff(ans2), c(ref2$u), tolerance=1e-8)
})

test_that("RRBLUP_D with fixed components reads the same markers that pullSnpGeno returns", {
  d = gsPop(nInd=60, nQtl=8, nSnp=12, nChr=3)
  y = d$pop@pheno[, 1L, drop=FALSE]
  G = pullSnpGeno(d$pop, simParam=d$SP)
  Ma = G - 1
  Md = domCode(G, 2)
  X = cbind(1, rowSums(Md))
  ref = solveRRBLUP_EM2(y, X, Ma, Md, Vu1=0.1, Vu2=0.05, Ve=1, tol=1e-6,
                        maxIter=10L, useEM=FALSE)
  ans = RRBLUP_D(d$pop, Vu=c(0.1, 0.05), Ve=1, simParam=d$SP)
  expect_equal(addEff(ans), c(ref$u[, 1L]), tolerance=1e-8)
})

test_that("RRBLUP_GCA and RRBLUP_SCA with fixed components read the parental haplotypes", {
  d = gsPop(nInd=60, nQtl=8, nSnp=12, nChr=3)
  y = d$pop@pheno[, 1L, drop=FALSE]
  X1 = matrix(1, nrow(y), 1L)
  M1 = 2*pullSnpHaplo(d$pop, haplo=1, simParam=d$SP) - 1
  M2 = 2*pullSnpHaplo(d$pop, haplo=2, simParam=d$SP) - 1
  ref = solveRRBLUP_EM2(y, X1, M1, M2, Vu1=0.1, Vu2=0.2, Ve=1, tol=1e-6,
                        maxIter=10L, useEM=FALSE)
  gca = RRBLUP_GCA(d$pop, Vu=c(0.1, 0.2), Ve=1,
                   simParam=d$SP)
  expect_equal(c(gca@gv[[1L]]@addEff), c(ref$u[, 1L]), tolerance=1e-8)
  expect_equal(c(gca@gv[[1L]]@addEffMale), c(ref$u[, 2L]), tolerance=1e-8)

  M3 = domCode(pullSnpGeno(d$pop, simParam=d$SP), 2)
  ref3 = solveRRBLUP_EM3(y, cbind(1, rowSums(M3)), M1, M2, M3, Vu1=0.1,
                         Vu2=0.2, Vu3=0.05, Ve=1, tol=1e-6, maxIter=10L,
                         useEM=FALSE)
  sca = RRBLUP_SCA(d$pop, Vu=c(0.1, 0.2, 0.05), Ve=1,
                   simParam=d$SP)
  expect_equal(c(sca@gv[[1L]]@addEff), c(ref3$u[, 1L]), tolerance=1e-8)
  expect_equal(c(sca@gv[[1L]]@addEffMale), c(ref3$u[, 2L]), tolerance=1e-8)
})

test_that("RRBLUP_GCA with fixed components splits a tetraploid's haplotypes in two", {
  set.seed(8401)
  founderPop = quickHaplo(nInd=60, nChr=2, segSites=20, ploidy=4L)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$restrSegSites(minQtlPerChr=10, minSnpPerChr=10, overlap=FALSE)
  SP$addTraitA(nQtlPerChr=10)
  SP$addSnpChip(nSnpPerChr=10)
  SP$setVarE(h2=0.5)
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)
  y = pop@pheno[, 1L, drop=FALSE]
  H = function(h) pullSnpHaplo(pop, haplo=h, simParam=SP)
  # Each half of the haplotypes is coded at half the ploidy level
  M1 = H(1) + H(2) - 1
  M2 = H(3) + H(4) - 1
  ref = solveRRBLUP_EM2(y, matrix(1, nrow(y), 1L), M1, M2, Vu1=0.1,
                        Vu2=0.2, Ve=1, tol=1e-6, maxIter=10L, useEM=FALSE)
  gca = RRBLUP_GCA(pop, Vu=c(0.1, 0.2), Ve=1, simParam=SP)
  expect_equal(c(gca@gv[[1L]]@addEff), c(ref$u[, 1L]), tolerance=1e-8)
  expect_equal(c(gca@gv[[1L]]@addEffMale), c(ref$u[, 2L]), tolerance=1e-8)
})

test_that("setEBV reads the parental haplotypes of GCA and SCA solutions", {
  skip_on_cran()
  d = gsPop(nInd=60, nQtl=8, nSnp=12, nChr=3)
  H1 = pullSnpHaplo(d$pop, haplo=1, simParam=d$SP)
  H2 = pullSnpHaplo(d$pop, haplo=2, simParam=d$SP)
  het = domCode(H1 + H2, 2)
  gcaSol = RRBLUP_GCA(d$pop, Vu=c(0.1, 0.2), Ve=1,
                      simParam=d$SP)
  gca = gcaSol@gv[[1L]]
  expected = gca@intercept + (2*H1 - 1)%*%gca@addEff +
    (2*H2 - 1)%*%gca@addEffMale
  expect_equal(c(ebv(setEBV(d$pop, gcaSol, simParam=d$SP))), c(expected),
               tolerance=1e-10)
  scaSol = RRBLUP_SCA(d$pop, Vu=c(0.1, 0.2, 0.05), Ve=1, simParam=d$SP)
  sca = scaSol@gv[[1L]]
  expected = sca@intercept + (2*H1 - 1)%*%sca@addEff +
    (2*H2 - 1)%*%sca@addEffMale + het%*%sca@domEff
  expect_equal(c(ebv(setEBV(d$pop, scaSol, simParam=d$SP))), c(expected),
               tolerance=1e-10)
})

test_that("setEBV with a targetPop uses the target's allele frequencies", {
  d = gsPop(nInd=60, nQtl=8, nSnp=12, nChr=3)
  target = d$pop[1:25]
  ans = RRBLUP(d$pop, simParam=d$SP)
  p = colMeans(pullSnpGeno(target, simParam=d$SP))/2
  a = addEff(ans)
  # With no dominance the breeding value is the additive effect centered
  # on the target population
  expected = (pullSnpGeno(d$pop, simParam=d$SP) - 1)%*%a - sum((2*p - 1)*a)
  out = setEBV(d$pop, ans, value="bv", targetPop=target, simParam=d$SP)
  expect_equal(c(ebv(out)), c(expected), tolerance=1e-10)
})

test_that("the dominance and hybrid models are not affected by the number of threads", {
  skip_on_cran()
  skip_if(getNumThreads() < 2L)
  set.seed(8501)
  founderPop = quickHaplo(nInd=60, nChr=2, segSites=20)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$restrSegSites(minQtlPerChr=10, minSnpPerChr=10, overlap=FALSE)
  SP$addTraitAD(nQtlPerChr=10, meanDD=0.5)
  SP$addSnpChip(nSnpPerChr=10)
  SP$setVarE(h2=0.5)
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)
  for(fit in list(RRBLUP_D, RRBLUP_GCA, RRBLUP_SCA)){
    one = fit(pop, nThreads=1L, simParam=SP)
    two = fit(pop, nThreads=2L, simParam=SP)
    expect_equal(addEff(one), addEff(two), tolerance=1e-8)
  }
})

# Fast mode is only taken when it is expected to be faster, which it never is
# on a population small enough to check against a dense solution. This takes
# it whatever it is expected to cost.
forceFast = function(expr){
  old = options(AlphaSimR.fastModeRoute="fast")
  on.exit(options(old))
  expr
}

# Multiple trait fast mode RRBLUP. With the variance components supplied it is
# the BLUP of the multiple trait model, which Henderson's equations give
# directly for a population small enough to write them out in R. With every
# record used for the variance components it is the same REML fit as multiple
# trait RRBLUP in its default mode, so the two have to agree.

# A population with three correlated traits
mvPop = function(nInd, nSnp, ploidy=2L, seed=8601){
  set.seed(seed)
  founderPop = quickHaplo(nInd=nInd, nChr=2, segSites=10+nSnp,
                          ploidy=ploidy)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$restrSegSites(minQtlPerChr=10, minSnpPerChr=nSnp, overlap=FALSE)
  SP$addTraitA(nQtlPerChr=10, mean=c(0,0,0), var=c(1,1,1),
               corA=matrix(c(1,0.6,-0.3,0.6,1,0.2,-0.3,0.2,1), nrow=3))
  SP$addSnpChip(nSnpPerChr=nSnp)
  SP$setVarE(h2=c(0.5,0.3,0.6))
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)
  pop@fixEff = rep_len(1:3, pop@nInd)
  return(list(pop=pop, SP=SP))
}

# Marker effects of every trait, one column per trait
addEffAll = function(ans){
  return(sapply(ans@gv, function(x) c(x@addEff)))
}

# Henderson's equations for the multiple trait model with known Vu and Ve,
# with the markers centered as fast mode RRBLUP centers them
denseMvBlup = function(pop, SP, Vu, Ve){
  P = pop@ploidy
  M = (pullSnpGeno(pop, simParam=SP) - P/2)*(2/P)
  Mc = scale(M, scale=FALSE)
  X = model.matrix(~factor(pop@fixEff))
  Y = pop@pheno
  nT = ncol(Y)
  W = cbind(kronecker(diag(nT), X), kronecker(diag(nT), Mc))
  Rinv = kronecker(solve(Ve), diag(nrow(Y)))
  C = crossprod(W, Rinv %*% W)
  fix = seq_len(ncol(X)*nT)
  C[-fix, -fix] = C[-fix, -fix] + kronecker(solve(Vu), diag(ncol(M)))
  sol = solve(C, crossprod(W, Rinv %*% c(Y)))
  return(matrix(sol[-fix], ncol(M), nT))
}

mvVu = matrix(c(0.02,0.01,-0.004,0.01,0.03,0.002,-0.004,0.002,0.025), 3)
mvVe = matrix(c(1,0.3,0.1,0.3,1.5,-0.2,0.1,-0.2,0.8), 3)

test_that("multiple trait fast mode RRBLUP is the BLUP of the multiple trait model", {
  d = mvPop(nInd=80, nSnp=15)
  ans = RRBLUP(d$pop, fastMode=TRUE, traits=1:3, Vu=mvVu, Ve=mvVe,
               simParam=d$SP)
  ref = denseMvBlup(d$pop, d$SP, mvVu, mvVe)
  # It iterates rather than solving, so it is held to a looser tolerance
  expect_equal(unname(addEffAll(ans)), ref, tolerance=1e-6)
  expect_equal(length(ans@bv), 3L)
  expect_equal(dimnames(ans@Vu), list(d$SP$traitNames, d$SP$traitNames))
})

test_that("multiple trait fast mode RRBLUP reads a tetraploid's markers", {
  # A ploidy other than two goes through the general loop of the products
  d = mvPop(nInd=60, nSnp=12, ploidy=4L)
  ans = RRBLUP(d$pop, fastMode=TRUE, traits=1:3, Vu=mvVu, Ve=mvVe,
               simParam=d$SP)
  ref = denseMvBlup(d$pop, d$SP, mvVu, mvVe)
  expect_equal(unname(addEffAll(ans)), ref, tolerance=1e-6)
})

test_that("multiple trait fast mode RRBLUP agrees with RRBLUP on every record", {
  skip_on_cran()
  # More records than markers decomposes the markers, and fewer decomposes
  # the record cross product built from the genotypes
  for(size in list(c(nInd=80, nSnp=15), c(nInd=40, nSnp=30))){
    d = mvPop(nInd=size[["nInd"]], nSnp=size[["nSnp"]])
    ans = forceFast(RRBLUP(d$pop, fastMode=TRUE, traits=1:3, nSubsample=0L,
                           simParam=d$SP))
    ref = RRBLUP(d$pop, traits=1:3, simParam=d$SP)
    info = paste("nInd", size[["nInd"]])
    expect_equal(unname(ans@Vu), unname(ref@Vu), tolerance=1e-8, info=info)
    expect_equal(unname(ans@Ve), unname(ref@Ve), tolerance=1e-8, info=info)
    expect_equal(addEffAll(ans), addEffAll(ref), tolerance=1e-6, info=info)
    # The genetic values are predicted from the uncentered markers, so the
    # intercepts have to agree as well as the effects
    expect_equal(c(ebv(setEBV(d$pop, ans, simParam=d$SP))),
                 c(ebv(setEBV(d$pop, ref, simParam=d$SP))),
                 tolerance=1e-6, info=info)
  }
})

test_that("fast mode RRBLUP predicts the same genetic values as RRBLUP", {
  # The intercept of the genetic values goes with the uncentered markers.
  # Taking it from the centered ones shifts every prediction by Mmean*u,
  # which is not zero in general, so this checks a trait where it is not.
  d = mvPop(nInd=80, nSnp=15)
  ans = forceFast(RRBLUP(d$pop, fastMode=TRUE, traits=1, nSubsample=0L,
                         simParam=d$SP))
  ref = RRBLUP(d$pop, traits=1, simParam=d$SP)
  expect_gt(abs(ref@bv[[1]]@intercept), 1e-8)
  expect_equal(ans@gv[[1]]@intercept, ref@gv[[1]]@intercept, tolerance=1e-6)
  expect_equal(c(ebv(setEBV(d$pop, ans, simParam=d$SP))),
               c(ebv(setEBV(d$pop, ref, simParam=d$SP))), tolerance=1e-6)
})

test_that("multiple trait fast mode RRBLUP leaves its inputs unchanged", {
  d = mvPop(nInd=40, nSnp=10)
  Vu = mvVu + 0
  Ve = mvVe + 0
  pheno0 = d$pop@pheno + 0
  RRBLUP(d$pop, fastMode=TRUE, traits=1:3, Vu=Vu, Ve=Ve, simParam=d$SP)
  expect_identical(Vu, mvVu)
  expect_identical(Ve, mvVe)
  expect_identical(d$pop@pheno, pheno0)
})

test_that("multiple trait fast mode RRBLUP rejects unusable variance components", {
  d = mvPop(nInd=40, nSnp=10)
  expect_error(RRBLUP(d$pop, fastMode=TRUE, traits=1:3, Vu=diag(2),
                      Ve=diag(3), simParam=d$SP), "one row per trait")
  expect_error(RRBLUP(d$pop, fastMode=TRUE, traits=1:3, Vu=mvVu,
                      Ve=diag(c(1,1,0)), simParam=d$SP), "positive definite")
  expect_error(RRBLUP(d$pop, fastMode=TRUE, traits=1:3, Vu=-mvVu, Ve=mvVe,
                      simParam=d$SP), "positive semidefinite")
  expect_error(RRBLUP(d$pop, fastMode=TRUE, traits=1:3, Vu=0*mvVu, Ve=mvVe,
                      simParam=d$SP), "not positive")
})

test_that("multiple trait fast mode RRBLUP rejects an indefinite Vu", {
  d = mvPop(nInd=40, nSnp=10)
  expect_error(RRBLUP(d$pop, fastMode=TRUE, traits=1:3,
                      Vu=diag(c(0.01,0.01,-0.005)), Ve=mvVe, simParam=d$SP),
               "positive semidefinite")
})

test_that("multiple trait fast mode RRBLUP checks Vu on its own scale", {
  # An error variance close to zero for the second trait would make the
  # first trait's negative variance look like rounding after whitening
  d = mvPop(nInd=40, nSnp=10)
  expect_error(RRBLUP(d$pop, fastMode=TRUE, traits=1:2,
                      Vu=diag(c(-0.01, 0.01)), Ve=diag(c(1, 1e-14)),
                      simParam=d$SP),
               "positive semidefinite")
})

test_that("multiple trait fast mode RRBLUP warns when maxIter is reached", {
  d = mvPop(nInd=80, nSnp=15)
  expect_warning(RRBLUP(d$pop, fastMode=TRUE, traits=1:3, Vu=mvVu, Ve=mvVe,
                        maxIter=1, simParam=d$SP), "without converging")
})

test_that("multiple trait fast mode RRBLUP with diagonal covariance matrices fits each trait on its own", {
  # With no covariance between the traits, the multiple trait BLUP of each
  # trait is its single trait BLUP. An error variance close to zero for the
  # first trait makes its variance ratio far larger than the second's, which
  # must not stop the second trait from being fitted or from converging.
  d = mvPop(nInd=80, nSnp=15)
  Vu = diag(c(0.02, 0.03))
  Ve = diag(c(1e-14, 1))
  ans = RRBLUP(d$pop, fastMode=TRUE, traits=1:2, Vu=Vu, Ve=Ve, maxIter=10000,
               simParam=d$SP)
  one = RRBLUP(d$pop, fastMode=TRUE, traits=2, Vu=Vu[2,2], Ve=Ve[2,2],
               maxIter=10000, simParam=d$SP)
  expect_gt(max(abs(addEff(one))), 0)
  expect_equal(addEffAll(ans)[, 2], addEff(one), tolerance=1e-6)
})

test_that("models name a response that is not one column per trait", {
  d = mvPop(nInd=60, nSnp=12)
  # Random values are a single column whatever traits are asked for
  for(fit in list(RRBLUP_D, RRBLUP_GCA, RRBLUP_SCA)){
    ans = fit(d$pop, traits=1:2, use="rand", simParam=d$SP)
    expect_equal(length(ans@gv), 1L)
    expect_equal(dim(ans@Ve), c(1L, 1L))
  }
  # A function of the traits returning two columns
  pair = function(x) x[, 1:2]
  ans = RRBLUP(d$pop, fastMode=TRUE, traits=pair, Vu=mvVu[1:2, 1:2],
               Ve=mvVe[1:2, 1:2], simParam=d$SP)
  expect_equal(dimnames(ans@Vu), list(c("Response1","Response2"),
                                      c("Response1","Response2")))
  expect_equal(ans@gv[[2]]@name, "est_GV_Response2")
  ans = RRBLUP_D(d$pop, traits=pair, simParam=d$SP)
  expect_equal(rownames(ans@Vu)[1:2], c("a:Response1", "a:Response2"))
})

test_that("models accept a function of the traits that returns a vector", {
  d = mvPop(nInd=60, nSnp=12)
  index = function(x) rowSums(x)
  ans = RRBLUP(d$pop, fastMode=TRUE, traits=index, Vu=0.01, Ve=1,
               simParam=d$SP)
  expect_equal(length(ans@gv), 1L)
  ans = RRBLUP(d$pop, fastMode=TRUE, traits=index, simParam=d$SP)
  expect_true(all(is.finite(addEff(ans))))
  for(fit in list(RRBLUP_D, RRBLUP_GCA, RRBLUP_SCA)){
    ans = fit(d$pop, traits=index, simParam=d$SP)
    expect_equal(length(ans@gv), 1L)
    expect_equal(dim(ans@Ve), c(1L, 1L))
  }
})

test_that("multiple trait fast mode RRBLUP needs every trait recorded", {
  expect_error(AlphaSimR:::checkCompleteRecords(matrix(c(1,NA,3,4), 2)),
               "every trait recorded")
  expect_null(AlphaSimR:::checkCompleteRecords(matrix(1:4, 2)))
})

test_that("multiple trait fast mode RRBLUP is not affected by the number of threads", {
  skip_on_cran()
  skip_if(getNumThreads() < 2L)
  d = mvPop(nInd=80, nSnp=15)
  one = forceFast(RRBLUP(d$pop, fastMode=TRUE, traits=1:3, nThreads=1L,
                         simParam=d$SP))
  two = forceFast(RRBLUP(d$pop, fastMode=TRUE, traits=1:3, nThreads=2L,
                         simParam=d$SP))
  expect_identical(addEffAll(one), addEffAll(two))
  expect_identical(one@Vu, two@Vu)
})

# Fast mode of the multikernel models. With the variance components supplied
# it is the BLUP of the model the direct form fits, which Henderson's
# equations give directly for a population small enough to write them out in
# R. With every record used for the variance components it is the same REML
# fit as the direct form, so the two have to agree.

mkModels = list(RRBLUP_D=RRBLUP_D, RRBLUP_GCA=RRBLUP_GCA,
                RRBLUP_SCA=RRBLUP_SCA)

# Each random term's codes, as the direct form of the model builds them
mkCodes = function(model, pop, SP){
  P = pop@ploidy
  x = pullSnpGeno(pop, simParam=SP)
  dom = x*(P-x)*(2/P)^2
  if(model=="RRBLUP_D"){
    return(lapply(list((x-P/2)*(2/P), dom), unname))
  }
  H = pullSnpHaplo(pop, simParam=SP)
  hap = rep(seq_len(P), pop@nInd)
  ind = rep(seq_len(pop@nInd), each=P)
  half = P/2
  code = function(take){
    (rowsum(H[take, , drop=FALSE], ind[take]) - half/2)*(2/half)
  }
  codes = list(code(hap<=half), code(hap>half))
  if(model=="RRBLUP_SCA"){
    codes[[3]] = dom
  }
  return(lapply(codes, unname))
}

# The fixed effects of the direct form: the sum to zero design of makeX, and
# for the models with dominance the summed dominance codes
mkFixed = function(model, pop, codes){
  f = factor(pop@fixEff)
  X = model.matrix(~f, contrasts.arg=list(f="contr.sum"))
  if(model!="RRBLUP_GCA"){
    X = cbind(X, rowSums(codes[[length(codes)]]))
  }
  return(unname(X))
}

# Henderson's equations for the multikernel model with known Vu and Ve
denseMkBlup = function(Y, X, Mlist, Vu, Ve){
  Y = as.matrix(Y)
  nT = ncol(Y)
  W = cbind(kronecker(diag(nT), X),
            do.call(cbind, lapply(Mlist, function(M) kronecker(diag(nT), M))))
  Rinv = kronecker(solve(Ve), diag(nrow(Y)))
  C = crossprod(W, Rinv %*% W)
  idx = list()
  at = ncol(X)*nT
  for(k in seq_along(Mlist)){
    idx[[k]] = at+seq_len(ncol(Mlist[[k]])*nT)
    at = max(idx[[k]])
    C[idx[[k]], idx[[k]]] = C[idx[[k]], idx[[k]]] +
      kronecker(solve(Vu[[k]]), diag(ncol(Mlist[[k]])))
  }
  sol = solve(C, crossprod(W, Rinv %*% c(Y)))
  beta = matrix(sol[seq_len(ncol(X)*nT)], ncol(X), nT)
  u = lapply(seq_along(Mlist),
             function(k) matrix(sol[idx[[k]]], ncol(Mlist[[k]]), nT))
  return(list(beta=beta, u=u))
}

# The effects a fit returns that come straight from the model, one matrix per
# term with a column per trait
mkEffects = function(model, ans){
  get = function(s, what){
    unname(sapply(slot(ans, s), function(x) slot(x, what)))
  }
  if(model=="RRBLUP_D"){
    return(list(get("gv", "addEff"), get("gv", "domEff")))
  }
  if(model=="RRBLUP_GCA"){
    return(list(get("female", "addEff"), get("male", "addEff")))
  }
  return(list(get("gv", "addEff"), get("gv", "addEffMale"),
              get("gv", "domEff")))
}

# What mkEffects should return for a dense solution. Dominance effects
# include the effect of the summed dominance codes, and are zero at markers
# where every individual has the same genotype.
mkExpected = function(model, ref, codes){
  eff = lapply(ref$u, function(u) unname(as.matrix(u)))
  if(model!="RRBLUP_GCA"){
    k = length(eff)
    meanD = ref$beta[nrow(ref$beta), ]
    if(model=="RRBLUP_SCA"){
      dose = codes[[1]]+codes[[2]]
    }else{
      dose = codes[[1]]
    }
    fixed = apply(dose, 2, function(x) all(x==x[1]))
    eff[[k]] = sweep(eff[[k]], 2, meanD, "+")
    eff[[k]][fixed, ] = 0
  }
  return(eff)
}

mkVu = function(model, nTraits){
  nTerm = if(model=="RRBLUP_SCA") 3 else 2
  scale = c(1, 0.5, 0.3)[seq_len(nTerm)]
  if(nTraits==1){
    return(0.02*scale)
  }
  return(lapply(scale, function(s) s*mvVu))
}

test_that("fast mode multikernel models are the BLUP of their models", {
  for(model in names(mkModels)){
    for(ploidy in c(2L, 4L)){
      if( (model=="RRBLUP_SCA") & (ploidy!=2L) ) next
      d = mvPop(nInd=80, nSnp=15, ploidy=ploidy)
      codes = mkCodes(model, d$pop, d$SP)
      X = mkFixed(model, d$pop, codes)
      for(traits in list(1L, 1:3)){
        nT = length(traits)
        Vu = mkVu(model, nT)
        Ve = if(nT==1) 1 else mvVe
        ans = mkModels[[model]](d$pop, traits=traits, fastMode=TRUE, Vu=Vu,
                                Ve=Ve, simParam=d$SP)
        VuList = if(nT==1) as.list(Vu) else Vu
        ref = denseMkBlup(d$pop@pheno[, traits, drop=FALSE], X, codes,
                          lapply(VuList, as.matrix), as.matrix(Ve))
        info = paste(model, "ploidy", ploidy, "traits", nT)
        eff = mkEffects(model, ans)
        expected = mkExpected(model, ref, codes)
        for(k in seq_along(eff)){
          # It iterates rather than solving, so it is held to a looser
          # tolerance
          expect_equal(matrix(eff[[k]], ncol=nT), expected[[k]],
                       tolerance=1e-6, info=paste(info, "term", k))
        }
        # The intercept goes with the uncentered codes, as in the direct form
        expect_equal(sapply(ans@gv, function(x) x@intercept), ref$beta[1, ],
                     tolerance=1e-6, info=info)
      }
    }
  }
})

test_that("fast mode multikernel models agree with direct mode on every record", {
  skip_on_cran()
  # More records than markers gives the REML fit the codes, and fewer gives
  # it a factor of their cross product
  for(size in list(c(nInd=80, nSnp=15), c(nInd=40, nSnp=30))){
    d = mvPop(nInd=size[["nInd"]], nSnp=size[["nSnp"]])
    for(model in names(mkModels)){
      for(traits in list(1L, 1:2)){
        info = paste(model, "nInd", size[["nInd"]], "traits", length(traits))
        ans = forceFast(mkModels[[model]](d$pop, traits=traits,
                                          fastMode=TRUE, nSubsample=0L,
                                          simParam=d$SP))
        ref = mkModels[[model]](d$pop, traits=traits, simParam=d$SP)
        expect_equal(ans@Vu, ref@Vu, tolerance=1e-6, info=info)
        expect_equal(ans@Ve, ref@Ve, tolerance=1e-6, info=info)
        if(model=="RRBLUP_D"){
          values = c("gv", "bv")
        }else{
          values = c("gv", "female", "male")
        }
        for(value in values){
          expect_equal(setEBV(d$pop, ans, value=value, onlyValues=TRUE,
                              simParam=d$SP),
                       setEBV(d$pop, ref, value=value, onlyValues=TRUE,
                              simParam=d$SP),
                       tolerance=1e-5, info=paste(info, value))
        }
      }
    }
  }
})

test_that("fast mode multikernel models take back their own variance components", {
  d = mvPop(nInd=60, nSnp=12)
  for(model in names(mkModels)){
    for(traits in list(1L, 1:2)){
      fit = forceFast(mkModels[[model]](d$pop, traits=traits,
                                        fastMode=TRUE, simParam=d$SP))
      again = mkModels[[model]](d$pop, traits=traits, fastMode=TRUE,
                                Vu=fit@Vu, Ve=fit@Ve, simParam=d$SP)
      info = paste(model, length(traits))
      expect_equal(again@Vu, fit@Vu, info=info)
      expect_equal(mkEffects(model, again), mkEffects(model, fit),
                   tolerance=1e-8, info=info)
    }
  }
})

test_that("the REML fit on a subset is the same from cross products as from codes", {
  # With at least as many markers as records in the subset, fast mode gives
  # the REML fit each term's cross product in place of its codes
  d = mvPop(nInd=80, nSnp=30)
  nSub = 50L
  for(model in names(mkModels)){
    nTerm = if(model=="RRBLUP_SCA") 3L else 2L
    codes = mkCodes(model, d$pop, d$SP)
    X = mkFixed(model, d$pop, codes)
    for(traits in list(1L, 1:2)){
      info = paste(model, length(traits))
      set.seed(1)
      ans = forceFast(mkModels[[model]](d$pop, traits=traits, fastMode=TRUE,
                                        nSubsample=nSub, simParam=d$SP))
      set.seed(1)
      s = AlphaSimR:::fastSubset(nSub, d$pop@nInd)
      ref = solveRRBLUPMK(pheno(d$pop)[s, traits, drop=FALSE],
                          X[s, , drop=FALSE],
                          lapply(codes, function(M) M[s, , drop=FALSE]))
      got = AlphaSimR:::fastMKVarComp(ans@Vu, ans@Ve, nTerm, length(traits))
      expect_equal(got$Vu, lapply(ref$Vu, function(v) unname(as.matrix(v))),
                   tolerance=1e-6, info=info)
      expect_equal(got$Ve, unname(as.matrix(ref$Ve)), tolerance=1e-6,
                   info=info)
    }
  }
})

test_that("the REML fit on a subset handles a low rank cross product", {
  # Many copies of two doubled haploid genotypes give the female and male
  # terms identical cross products of rank two. Their factors have two
  # columns, so the fit belongs on the coefficient side, and its average
  # information is singular, so it leans on EM, which divides each update by
  # a bound for the rank. The fit has to match the one given a factor of each
  # cross product with as many columns as its rank, as an eigendecomposition
  # gives. With identical cross products only the sum of the two variances
  # is identified, so that is what is compared.
  set.seed(8601)
  founderPop = quickHaplo(nInd=20, nChr=2, segSites=40)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$restrSegSites(minQtlPerChr=10, minSnpPerChr=30, overlap=FALSE)
  SP$addTraitA(nQtlPerChr=10, mean=c(0,0), var=c(1,1),
               corA=matrix(c(1,0.5,0.5,1), nrow=2))
  SP$addSnpChip(nSnpPerChr=30)
  SP$setVarE(h2=c(0.5,0.5))
  dh = makeDH(newPop(founderPop, simParam=SP)[1:2], nDH=1, simParam=SP)
  pop = setPheno(dh[rep(1:2, 40)], simParam=SP)
  nSub = 50L
  for(model in c("RRBLUP_GCA", "RRBLUP_SCA")){
    nTerm = if(model=="RRBLUP_SCA") 3L else 2L
    set.seed(1)
    s = AlphaSimR:::fastSubset(nSub, pop@nInd)
    codes = mkCodes(model, pop, SP)
    # One fixed effect level, so the intercept, and for RRBLUP_SCA the
    # summed dominance codes, which mkFixed would add
    X = matrix(1, nSub, 1)
    if(model=="RRBLUP_SCA"){
      X = cbind(X, rowSums(codes[[3]][s, , drop=FALSE]))
    }
    factors = lapply(codes, function(M){
      e = eigen(tcrossprod(M[s, , drop=FALSE]), symmetric=TRUE)
      keep = e$values > nSub*.Machine$double.eps*max(e$values)
      e$vectors[, keep, drop=FALSE] %*% diag(sqrt(e$values[keep]), sum(keep))
    })
    for(traits in list(1L, 1:2)){
      info = paste(model, length(traits))
      set.seed(1)
      ans = forceFast(mkModels[[model]](pop, traits=traits, fastMode=TRUE,
                                        nSubsample=nSub, simParam=SP))
      ref = solveRRBLUPMK(pheno(pop)[s, traits, drop=FALSE], X, factors,
                          maxIter=100)
      got = AlphaSimR:::fastMKVarComp(ans@Vu, ans@Ve, nTerm, length(traits))
      # For one trait ref$Vu is a matrix with a row per term, which lapply
      # takes element by element
      refVu = lapply(ref$Vu, function(v) unname(as.matrix(v)))
      expect_equal(got$Vu[[1]]+got$Vu[[2]], refVu[[1]]+refVu[[2]],
                   tolerance=1e-6, info=info)
      expect_equal(got$Ve, unname(as.matrix(ref$Ve)), tolerance=1e-6,
                   info=info)
    }
  }
})

test_that("each model records its name in the fit", {
  d = mvPop(nInd=40, nSnp=10)
  expect_identical(RRBLUP(d$pop, simParam=d$SP)@model, "RRBLUP")
  expect_identical(forceFast(RRBLUP(d$pop, fastMode=TRUE,
                                    simParam=d$SP))@model, "RRBLUP")
  # The numbered models fit the same models, so they record the same names.
  # They are given their variance components, so that EM is not run.
  expect_identical(quietDeprecated(RRBLUP2(d$pop, Vu=0.01, Ve=1, useEM=FALSE,
                                           simParam=d$SP))@model, "RRBLUP")
  numbered = quietDeprecated(list(
    RRBLUP_D=RRBLUP_D2(d$pop, Va=0.01, Vd=0.01, Ve=1, useEM=FALSE,
                       simParam=d$SP),
    RRBLUP_GCA=RRBLUP_GCA2(d$pop, VuF=0.01, VuM=0.01, Ve=1, useEM=FALSE,
                           simParam=d$SP),
    RRBLUP_SCA=RRBLUP_SCA2(d$pop, VuF=0.01, VuM=0.01, VuD=0.01, Ve=1,
                           useEM=FALSE, simParam=d$SP)))
  for(model in names(mkModels)){
    expect_identical(mkModels[[model]](d$pop, simParam=d$SP)@model, model)
    expect_identical(forceFast(mkModels[[model]](d$pop, fastMode=TRUE,
                                                 simParam=d$SP))@model,
                     model)
    expect_identical(numbered[[model]]@model, model, info=model)
  }
  expect_identical(new("RRsol")@model, "")
})

test_that("fast mode takes the variance components of an earlier fit", {
  d = mvPop(nInd=60, nSnp=12)
  for(traits in list(1L, 1:2)){
    fit = RRBLUP(d$pop, traits=traits, simParam=d$SP)
    expect_identical(RRBLUP(d$pop, traits=traits, fastMode=TRUE,
                            varComp=fit, simParam=d$SP),
                     RRBLUP(d$pop, traits=traits, fastMode=TRUE,
                            Vu=fit@Vu, Ve=fit@Ve, simParam=d$SP),
                     info=length(traits))
  }
  # A fit of the numbered model is a fit of the same model
  fit2 = quietDeprecated(RRBLUP2(d$pop, Vu=0.01, Ve=1, useEM=FALSE,
                                 simParam=d$SP))
  expect_identical(RRBLUP(d$pop, fastMode=TRUE, varComp=fit2, simParam=d$SP),
                   RRBLUP(d$pop, fastMode=TRUE, Vu=fit2@Vu, Ve=fit2@Ve,
                          simParam=d$SP))
  for(model in names(mkModels)){
    for(traits in list(1L, 1:2)){
      fit = mkModels[[model]](d$pop, traits=traits, simParam=d$SP)
      expect_identical(mkModels[[model]](d$pop, traits=traits,
                                         fastMode=TRUE, varComp=fit,
                                         simParam=d$SP),
                       mkModels[[model]](d$pop, traits=traits,
                                         fastMode=TRUE, Vu=fit@Vu,
                                         Ve=fit@Ve, simParam=d$SP),
                       info=paste(model, length(traits)))
    }
  }
})

test_that("an earlier fit is refused when it does not match the model", {
  d = mvPop(nInd=40, nSnp=10)
  fitA = RRBLUP(d$pop, simParam=d$SP)
  fitD = RRBLUP_D(d$pop, simParam=d$SP)
  expect_error(RRBLUP(d$pop, fastMode=TRUE, varComp=fitA, Vu=0.01,
                      simParam=d$SP), "not both")
  expect_error(RRBLUP(d$pop, varComp=fitA, Vu=0.01, simParam=d$SP),
               "not both")
  expect_error(RRBLUP_D(d$pop, varComp=fitA, simParam=d$SP),
               "fitted with RRBLUP, not RRBLUP_D")
  expect_error(RRBLUP(d$pop, fastMode=TRUE, varComp=fitD, simParam=d$SP),
               "fitted with RRBLUP_D, not RRBLUP")
  # Two models with two terms each, which the values alone cannot tell apart
  expect_error(RRBLUP_GCA(d$pop, fastMode=TRUE, varComp=fitD,
                          simParam=d$SP),
               "fitted with RRBLUP_D, not RRBLUP_GCA")
  expect_error(RRBLUP(d$pop, traits=1:2, fastMode=TRUE, varComp=fitA,
                      simParam=d$SP), "fitted to 1 trait")
  # The first trait's QTL, as many as the markers but at other loci
  expect_error(RRBLUP(d$pop, fastMode=TRUE, varComp=fitA, useQtl=TRUE,
                      snpChip=1, simParam=d$SP), "different markers")
  expect_error(RRBLUP(d$pop, fastMode=TRUE, varComp=list(), simParam=d$SP),
               "class RRsol")
  noModel = fitA
  noModel@model = ""
  expect_error(RRBLUP(d$pop, fastMode=TRUE, varComp=noModel, simParam=d$SP),
               "does not record")
  # A fit saved before the slot existed has no model attribute at all
  old = fitA
  attr(old, "model") = NULL
  expect_error(RRBLUP(d$pop, fastMode=TRUE, varComp=old, simParam=d$SP),
               "does not record")
})

test_that("a term with no variance has no effects in fast mode", {
  d = mvPop(nInd=60, nSnp=12)
  ans = RRBLUP_D(d$pop, fastMode=TRUE, Vu=c(0.02, 0), Ve=1, simParam=d$SP)
  dom = ans@gv[[1]]@domEff
  # Every dominance effect is the effect of the summed dominance codes, or
  # zero where the marker is fixed
  expect_lte(length(unique(dom[dom!=0])), 1L)
  expect_gt(max(abs(ans@gv[[1]]@addEff)), 0)
})

test_that("fast mode multikernel models check their arguments", {
  d = mvPop(nInd=40, nSnp=10)
  for(model in names(mkModels)){
    fit = mkModels[[model]]
    nTerm = if(model=="RRBLUP_SCA") 3 else 2
    # Outside fast mode, estimating Ve alone would ignore Vu without notice
    expect_error(fit(d$pop, Vu=rep(0.01, nTerm), simParam=d$SP),
                 "both Vu and Ve", info=model)
    expect_error(fit(d$pop, fastMode=TRUE, Vu=0.01, Ve=1, simParam=d$SP),
                 "one value per random term", info=model)
    expect_error(fit(d$pop, traits=1:2, fastMode=TRUE,
                     Vu=matrix(0.01, 2*nTerm, 2*nTerm),
                     Ve=diag(2), simParam=d$SP),
                 "block diagonal", info=model)
    Vu = rep(list(diag(2)), nTerm)
    Vu[[1]] = diag(c(0.01, -0.01))
    expect_error(fit(d$pop, traits=1:2, fastMode=TRUE, Vu=Vu, Ve=diag(2),
                     simParam=d$SP),
                 "positive semidefinite", info=model)
    expect_error(fit(d$pop, traits=1:2, fastMode=TRUE,
                     Vu=rep(list(diag(2)), nTerm), Ve=diag(c(1, 0)),
                     simParam=d$SP),
                 "positive definite", info=model)
    expect_warning(fit(d$pop, fastMode=TRUE, Vu=rep(0.01, nTerm), Ve=1,
                       maxIter=1L, simParam=d$SP),
                   "without converging", info=model)
  }
  d4 = mvPop(nInd=40, nSnp=10, ploidy=4L)
  expect_error(forceFast(RRBLUP_SCA(d4$pop, fastMode=TRUE, simParam=d4$SP)),
               "only diploids")
})

test_that("fast mode multikernel models are not affected by the number of threads", {
  skip_on_cran()
  skip_if(getNumThreads() < 2L)
  d = mvPop(nInd=80, nSnp=15)
  for(model in names(mkModels)){
    for(traits in list(1L, 1:3)){
      one = forceFast(mkModels[[model]](d$pop, traits=traits,
                                        fastMode=TRUE, nThreads=1L,
                                        simParam=d$SP))
      two = forceFast(mkModels[[model]](d$pop, traits=traits,
                                        fastMode=TRUE, nThreads=2L,
                                        simParam=d$SP))
      info = paste(model, length(traits))
      expect_identical(mkEffects(model, one), mkEffects(model, two),
                       info=info)
      expect_identical(one@Vu, two@Vu, info=info)
    }
  }
})

test_that("RRBLUPMemUse counts the fast mode of the multikernel models", {
  for(m in c("RRBLUP_D","RRBLUP_GCA","RRBLUP_SCA")){
    small = RRBLUPMemUse(nInd=500, nMarker=500, model=m, fastMode=TRUE)
    large = RRBLUPMemUse(nInd=20000, nMarker=5000, model=m, fastMode=TRUE)
    expect_true(is.finite(small), info=m)
    expect_gt(large, small)
    # Far less than the direct form once the records are many
    expect_lt(large, RRBLUPMemUse(nInd=20000, nMarker=5000, model=m),
              label=m)
    expect_gt(RRBLUPMemUse(nInd=20000, nMarker=5000, model=m, nTraits=3,
                           fastMode=TRUE), large)
  }
})

test_that("fast mode is taken only when it is expected to be faster", {
  route = AlphaSimR:::useFastRoute
  # Supplied variance components need the fast route
  expect_true(route("RRBLUP_D", 0.1, 1, 100, 50, 1, 1, 5000))
  # A subset of every record would repeat the direct fit
  expect_false(route("RRBLUP_D", NULL, NULL, 4000, 50000, 1, 1, 5000))
  expect_false(route("RRBLUP_D", NULL, NULL, 50000, 50000, 1, 1, 0))
  # Many records and many markers favor the fast route
  for(m in c("RRBLUP", "RRBLUP_D", "RRBLUP_GCA", "RRBLUP_SCA")){
    expect_true(route(m, NULL, NULL, 50000, 20000, 1, 1, 2000), info=m)
  }
  # Few markers make the direct fit of RRBLUP cheap
  expect_false(route("RRBLUP", NULL, NULL, 20000, 500, 1, 1, 2000))
})

test_that("fast mode that is not expected to be faster fits the direct form", {
  d = mvPop(nInd=60, nSnp=12)
  for(model in c("RRBLUP", names(mkModels))){
    fit = get(model)
    expect_identical(fit(d$pop, fastMode=TRUE, simParam=d$SP),
                     fit(d$pop, simParam=d$SP), info=model)
  }
})

test_that("fast mode multikernel models fit the fixed effects alone when every variance is zero", {
  d = mvPop(nInd=60, nSnp=12)
  y = d$pop@pheno[, 1]
  for(model in names(mkModels)){
    nTerm = if(model=="RRBLUP_SCA") 3 else 2
    codes = mkCodes(model, d$pop, d$SP)
    X = mkFixed(model, d$pop, codes)
    b = qr.solve(X, y)
    for(traits in list(1L, 1:2)){
      nT = length(traits)
      if(nT==1){
        Vu = rep(0, nTerm)
        Ve = 1
      }else{
        Vu = rep(list(matrix(0, 2, 2)), nTerm)
        Ve = diag(2)
      }
      ans = mkModels[[model]](d$pop, traits=traits, fastMode=TRUE, Vu=Vu,
                              Ve=Ve, simParam=d$SP)
      info = paste(model, nT)
      eff = mkEffects(model, ans)
      ref = list(beta=matrix(b, ncol=1),
                 u=lapply(codes, function(M) matrix(0, ncol(M), 1)))
      expected = mkExpected(model, ref, codes)
      for(k in seq_along(eff)){
        expect_equal(matrix(eff[[k]], ncol=nT)[, 1], expected[[k]][, 1],
                     tolerance=1e-8, info=paste(info, "term", k))
      }
      expect_equal(ans@gv[[1]]@intercept, b[1], tolerance=1e-8, info=info)
    }
  }
})

# At variance ratios near the precision of a double, a fit may warn that
# rounding limits its accuracy or that its iteration stalled. Those warnings
# describe the whole solution, while the tests below check the effects they
# can verify, so only those two warnings are let through quietly; any other,
# such as reaching maxIter, still fails expect_warning(..., NA) or is
# reported.
quietRounding = function(expr){
  withCallingHandlers(expr, warning=function(w){
    if(grepl("rounding limits|stalled", conditionMessage(w))){
      invokeRestart("muffleWarning")
    }
  })
}

test_that("fast mode multikernel models converge for every trait whatever their error variances", {
  # With no covariance between the traits, each trait's BLUP is its single
  # trait BLUP. An error variance close to zero for the first trait makes its
  # right hand sides far larger than the second's, which must not stop the
  # second trait from converging.
  d = mvPop(nInd=80, nSnp=15)
  for(model in names(mkModels)){
    nTerm = if(model=="RRBLUP_SCA") 3 else 2
    scale = c(1, 0.5, 0.3)[seq_len(nTerm)]
    Vu = lapply(scale, function(s) diag(c(0.02, 0.03))*s)
    Ve = diag(c(1e-14, 1))
    ans = quietRounding(mkModels[[model]](d$pop, traits=1:2, fastMode=TRUE,
                                          Vu=Vu, Ve=Ve, maxIter=10000L,
                                          simParam=d$SP))
    one = mkModels[[model]](d$pop, traits=2, fastMode=TRUE,
                            Vu=0.03*scale, Ve=1, maxIter=10000L,
                            simParam=d$SP)
    two = mkEffects(model, ans)
    single = mkEffects(model, one)
    for(k in seq_along(two)){
      expect_gt(max(abs(single[[k]])), 0)
      expect_equal(two[[k]][, 2], c(single[[k]]), tolerance=1e-6,
                   info=paste(model, "term", k))
    }
  }
})

test_that("the RR-BLUP models keep the positions of their arguments", {
  # The fast mode arguments come after the dots, so calls that give simParam
  # and nThreads by position are unchanged
  d = mvPop(nInd=60, nSnp=12)
  expect_identical(RRBLUP(d$pop, 1, "pheno", 1, FALSE, 1000L, d$SP, 1L),
                   RRBLUP(d$pop, simParam=d$SP, nThreads=1L))
  for(model in names(mkModels)){
    fit = mkModels[[model]]
    expect_identical(fit(d$pop, 1, "pheno", 1, FALSE, 100L, d$SP, 1L),
                     fit(d$pop, simParam=d$SP, nThreads=1L), info=model)
  }
})

test_that("fast mode converges when a coupled term has no right hand side", {
  # A response orthogonal to the paternal codes and the fixed effects gives
  # that term a right hand side of zero, while the maternal term, which is
  # coupled to it through the records, has one. The paternal term's
  # tolerance is then a fraction of nothing, which rounding cannot reach, and
  # the fit must still stop without a warning and give the BLUP.
  d = mvPop(nInd=40, nSnp=4)
  codes = mkCodes("RRBLUP_GCA", d$pop, d$SP)
  X = mkFixed("RRBLUP_GCA", d$pop, codes)
  set.seed(31)
  z = rnorm(d$pop@nInd)
  y = qr.resid(qr(cbind(X, codes[[2]])), z)
  expect_lt(max(abs(crossprod(codes[[2]], y))), 1e-12)
  d$pop@pheno[, 1] = y+10
  Vu = c(0.02, 0.02)
  expect_warning(ans <- RRBLUP_GCA(d$pop, fastMode=TRUE, Vu=Vu, Ve=1,
                                   simParam=d$SP), NA)
  ref = denseMkBlup(d$pop@pheno[, 1, drop=FALSE], X, codes,
                    lapply(Vu, as.matrix), matrix(1))
  eff = mkEffects("RRBLUP_GCA", ans)
  expected = mkExpected("RRBLUP_GCA", ref, codes)
  for(k in seq_along(eff)){
    expect_equal(matrix(eff[[k]], ncol=1), expected[[k]], tolerance=1e-6,
                 info=paste("term", k))
  }
})

test_that("fast mode handles terms that coincide at every marker", {
  # In doubled haploids the maternal and paternal codes are the same, so the
  # two GCA terms coincide at every marker, and with a variance ratio this
  # large the identity in each block of the preconditioner is lost to
  # rounding. The fit must still succeed. With equal variances the two terms
  # share the effects equally, and together they are RRBLUP on the same codes
  # with twice the variance.
  d = mvPop(nInd=100, nSnp=12)
  dh = setPheno(makeDH(d$pop, simParam=d$SP), simParam=d$SP)
  ans = quietRounding(RRBLUP_GCA(dh, fastMode=TRUE, Vu=c(0.02, 0.02),
                                 Ve=1e-16, maxIter=10000L, simParam=d$SP))
  female = ans@female[[1]]@addEff
  male = ans@male[[1]]@addEff
  expect_true(all(is.finite(female)))
  expect_gt(max(abs(female)), 0)
  expect_equal(female, male, tolerance=1e-6)
  ref = RRBLUP(dh, fastMode=TRUE, Vu=0.04, Ve=1e-16, maxIter=10000L,
               simParam=d$SP)
  expect_equal(2*female, ref@bv[[1]]@addEff, tolerance=1e-6)
})

test_that("coinciding terms share their effects in proportion to their variances", {
  # Effects along the difference of two coinciding terms are unseen by the
  # records and held at zero by the prior alone, which a large variance
  # ratio dwarfs, so they are removed exactly. The terms then split each
  # effect of RRBLUP with their summed variance in proportion to their own.
  d = mvPop(nInd=100, nSnp=12)
  dh = setPheno(makeDH(d$pop, simParam=d$SP), simParam=d$SP)
  # Down to 1e-12, which every platform resolves; at 1e-16 whether the
  # iteration converges depends on the platform's rounding
  for(Ve in c(1, 1e-8, 1e-12)){
    ans = quietRounding(RRBLUP_GCA(dh, fastMode=TRUE, Vu=c(0.02, 0.03),
                                   Ve=Ve, maxIter=10000L, simParam=d$SP))
    ref = RRBLUP(dh, fastMode=TRUE, Vu=0.05, Ve=Ve, maxIter=10000L,
                 simParam=d$SP)
    u = ref@bv[[1]]@addEff
    expect_equal(ans@female[[1]]@addEff, 0.4*u, tolerance=1e-6,
                 info=paste("Ve", Ve))
    expect_equal(ans@male[[1]]@addEff, 0.6*u, tolerance=1e-6,
                 info=paste("Ve", Ve))
  }
  # Two traits, one with an error variance far below the other's. At 1e-12
  # the split drifts by about 1e-4 without removing the unseen directions,
  # while the solve stays well within double precision on every platform; at
  # 1e-16 whether it converges depends on the platform's rounding.
  V = diag(c(0.02, 0.03))
  Ve = diag(c(1e-12, 1))
  expect_warning(ans <- quietRounding(
    RRBLUP_GCA(dh, traits=1:2, fastMode=TRUE, Vu=list(V, 1.5*V), Ve=Ve,
               maxIter=10000L, simParam=d$SP)), NA)
  ref = RRBLUP(dh, traits=1:2, fastMode=TRUE, Vu=2.5*V, Ve=Ve,
               maxIter=10000L, simParam=d$SP)
  u = sapply(ref@bv, function(x) x@addEff)
  expect_equal(sapply(ans@female, function(x) x@addEff), 0.4*u,
               tolerance=1e-6)
  expect_equal(sapply(ans@male, function(x) x@addEff), 0.6*u,
               tolerance=1e-6)
  # Singular covariance matrices, with variance in the first trait alone.
  # The unseen directions are then confined to the trait directions both
  # terms carry.
  for(Ve in list(diag(2), diag(c(1e-16, 1)))){
    ans = quietRounding(
      RRBLUP_GCA(dh, traits=1:2, fastMode=TRUE,
                 Vu=list(diag(c(0.02, 0)), diag(c(0.03, 0))), Ve=Ve,
                 maxIter=10000L, simParam=d$SP))
    ref = RRBLUP(dh, traits=1:2, fastMode=TRUE, Vu=diag(c(0.05, 0)), Ve=Ve,
                 maxIter=10000L, simParam=d$SP)
    u = sapply(ref@bv, function(x) x@addEff)
    info = paste("Ve", Ve[1,1])
    expect_equal(sapply(ans@female, function(x) x@addEff), 0.4*u,
                 tolerance=1e-6, info=info)
    expect_equal(sapply(ans@male, function(x) x@addEff), 0.6*u,
                 tolerance=1e-6, info=info)
  }
})

# A population of diploids made of groups of identical individuals, each
# group given by its maternal and paternal haplotypes as strings of alleles,
# with no traits, so that its response can be set by hand
haploPop = function(groups, nPerGroup){
  bits = function(h) as.integer(strsplit(h, "")[[1]])
  H = do.call(rbind, lapply(groups, function(g){
    rbind(bits(g[1]), bits(g[2]))[rep(1:2, nPerGroup), , drop=FALSE]
  }))
  nLoci = ncol(H)
  genMap = data.frame(markerName=paste0("m", seq_len(nLoci)), chromosome=1,
                      position=seq(0, 0.5, length.out=nLoci))
  colnames(H) = genMap$markerName
  founderPop = importHaplo(H, genMap, ploidy=2L)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addSnpChip(nSnpPerChr=nLoci)
  pop = newPop(founderPop, simParam=SP)
  return(list(pop=pop, SP=SP))
}

# The BLUP of a multikernel model from the variance matrix of the records,
# which needs no inverse of Vu, so Vu may be singular
recordMkBlup = function(Y, X, Mlist, Vu, Ve){
  Y = as.matrix(Y)
  n = nrow(Y)
  V = kronecker(Ve, diag(n))
  for(k in seq_along(Mlist)){
    V = V + kronecker(Vu[[k]], tcrossprod(Mlist[[k]]))
  }
  XX = kronecker(diag(ncol(Y)), X)
  Vi = solve(V)
  b = solve(crossprod(XX, Vi %*% XX), crossprod(XX, Vi %*% c(Y)))
  r = Vi %*% (c(Y)-XX %*% b)
  u = lapply(seq_along(Mlist), function(k){
    matrix(kronecker(Vu[[k]], t(Mlist[[k]])) %*% r, ncol(Mlist[[k]]),
           ncol(Y))
  })
  return(list(beta=matrix(b, ncol(X)), u=u))
}

test_that("fast mode gives the BLUP when markers depend on each other", {
  # Every individual of a group is the same, and the maternal and paternal
  # codes of the two markers are linearly dependent across the groups, so
  # some effects cancel across markers. Only the prior holds them, against
  # data many orders of magnitude stronger, so a residual small next to the
  # right hand sides would leave them far off.
  d = haploPop(list(c("00", "11"), c("10", "10"), c("01", "00")), 32)
  y = rep(c(0, 1, 3), each=32)
  d$pop@pheno = matrix(y, ncol=1)
  ans = RRBLUP_GCA(d$pop, fastMode=TRUE, Vu=c(0.02, 0.03), Ve=1e-8,
                   maxIter=10000L, simParam=d$SP)
  codes = mkCodes("RRBLUP_GCA", d$pop, d$SP)
  X = matrix(1, length(y), 1)
  ref = denseMkBlup(matrix(y), X, codes, list(matrix(0.02), matrix(0.03)),
                    matrix(1e-8))
  expect_equal(ans@female[[1]]@addEff, c(ref$u[[1]]), tolerance=1e-6)
  expect_equal(ans@male[[1]]@addEff, c(ref$u[[2]]), tolerance=1e-6)
  # With a variance ratio this much larger, rounding alone moves the effects
  # that cancel across the markers by more than the tolerance, and the fit
  # must say so rather than report convergence
  expect_warning(RRBLUP_GCA(d$pop, fastMode=TRUE, Vu=c(0.02, 0.03),
                            Ve=1e-12, maxIter=10000L, simParam=d$SP),
                 "rounding limits the accuracy")
})

test_that("fast mode gives the BLUP when several null directions share a marker", {
  # At the first marker the three SCA codes are proportional once the
  # intercept and the summed dominance codes are absorbed, so their cross
  # products have two null vectors. With each term's variance confined to
  # different trait directions, no single null vector is carried by every
  # term it involves, but combinations of them are, and those are removed
  # together. The reference is the BLUP from the variance matrix of the
  # records, which takes the singular covariance matrices as they are.
  d = haploPop(list(c("00", "00"), c("10", "10"), c("11", "01")), 32)
  set.seed(77)
  Y = cbind(rep(c(0, 1, 3), each=32), rep(c(2, -1, 1), each=32))+
    matrix(rnorm(2*96, sd=1e-3), ncol=2)
  d$pop@pheno = Y
  Vu = list(diag(c(0.02, 0)), diag(c(0, 0.03)), matrix(0.04, 2, 2))
  Ve = diag(c(1e-6, 1e-6))
  ans = RRBLUP_SCA(d$pop, traits=1:2, fastMode=TRUE, Vu=Vu, Ve=Ve,
                   maxIter=10000L, simParam=d$SP)
  codes = mkCodes("RRBLUP_SCA", d$pop, d$SP)
  X = cbind(1, rowSums(codes[[3]]))
  ref = recordMkBlup(Y, X, codes, Vu, Ve)
  expect_equal(sapply(ans@gv, function(x) x@addEff), ref$u[[1]],
               tolerance=1e-6)
  expect_equal(sapply(ans@gv, function(x) x@addEffMale), ref$u[[2]],
               tolerance=1e-6)
})

test_that("fast mode separates terms whose covariance matrices are nearly aligned", {
  # In doubled haploids the two GCA terms have the same code at the one
  # marker, and each term's variance lies along one trait direction, (1,0)
  # and (1,d), so the records tell them apart only through differences of
  # order d. Next to the variance ratio those are below the precision of the
  # preconditioner's blocks, which must not stop the iteration from
  # searching along them. With u1 = a*(1,0) and u2 = b*(1,d), each of prior
  # variance 1/64, and 100 records of code -1 or 1, the BLUP solves
  #   (100/d^2)*(a+b)+64*a = 0
  #   (100/d^2)*(a+b)-100*(1-b)+64*b = 0
  # whose solution is below.
  d = 2^-27
  dd = haploPop(list(c("0", "0"), c("1", "1")), 50)
  x = rep(c(-1, 1), each=50)
  dd$pop@pheno = cbind(0, d*x)
  Vu = list(diag(c(1/64, 0)), tcrossprod(c(1, d))/64)
  ans = quietRounding(RRBLUP_GCA(dd$pop, traits=1:2, fastMode=TRUE, Vu=Vu,
                                 Ve=diag(2)*d^2, maxIter=10000L,
                                 simParam=dd$SP))
  k = 0.64*d^2
  b = 100*(1+k)/(164*(1+k)+64)
  a = -b/(1+k)
  female = sapply(ans@female, function(x) x@addEff)
  male = sapply(ans@male, function(x) x@addEff)
  expect_equal(female[1], a, tolerance=1e-6)
  expect_equal(male[1], b, tolerance=1e-6)
  expect_lt(abs(female[2]), 1e-12)
  expect_equal(male[2], b*d, tolerance=1e-6)
})

test_that("fast mode multikernel models are not affected by the mean of the response", {
  # The intercept absorbs any shift of the response, so the marker effects
  # must not change. A large shift keeps the response's variation in its low
  # order bits, which must not be rounded away before the mean is taken off.
  # The responses are rounded to multiples of 2^-20, so that shifted by 2^30
  # they are still exact and the two fits are given the same information.
  d = mvPop(nInd=60, nSnp=12)
  d$pop@pheno = round(d$pop@pheno*2^20)/2^20
  shifted = d$pop
  shifted@pheno = shifted@pheno+2^30
  expect_identical(shifted@pheno-2^30, d$pop@pheno)
  for(model in names(mkModels)){
    nTerm = if(model=="RRBLUP_SCA") 3 else 2
    for(traits in list(1L, 1:2)){
      if(length(traits)==1){
        Vu = 0.02*c(1, 0.5, 0.3)[seq_len(nTerm)]
        Ve = 0.3
      }else{
        Vu = lapply(c(1, 0.5, 0.3)[seq_len(nTerm)],
                    function(s) s*mvVu[1:2, 1:2])
        Ve = mvVe[1:2, 1:2]
      }
      ref = mkModels[[model]](d$pop, traits=traits, fastMode=TRUE, Vu=Vu,
                              Ve=Ve, simParam=d$SP)
      ans = mkModels[[model]](shifted, traits=traits, fastMode=TRUE, Vu=Vu,
                              Ve=Ve, simParam=d$SP)
      expect_equal(mkEffects(model, ans), mkEffects(model, ref),
                   tolerance=1e-7, info=paste(model, length(traits)))
    }
  }
})

test_that("fast mode keeps each fixed effect level's variation whatever the levels' means", {
  # Two fixed effect levels of 32 individuals, with maternal codes z of -1 or
  # 1 balanced within each. The first level's responses are all 2^60 and the
  # second's are z, so the level effects absorb the first level entirely and
  # the BLUP of the maternal effect is z'y/(z'z+Ve/Vu) = 32/(64+50) = 16/57.
  # Giving the second level -z instead reverses it. Taking one offset off
  # every level would round the second level's variation away next to 2^60.
  d = haploPop(list(c("0", "0"), c("1", "0"), c("0", "0"), c("1", "0")), 16)
  d$pop@fixEff = rep(1:2, each=32)
  z = rep(c(-1, 1, -1, 1), each=16)
  for(sgn in c(1, -1)){
    d$pop@pheno = matrix(c(rep(2^60, 32), sgn*z[33:64]), ncol=1)
    ans = RRBLUP_GCA(d$pop, fastMode=TRUE, Vu=c(0.02, 0), Ve=1,
                     simParam=d$SP)
    expect_equal(ans@female[[1]]@addEff, sgn*16/57, tolerance=1e-8,
                 info=paste("sign", sgn))
  }
})

test_that("DEPRECATED numbered models warn and match their replacements", {
  # Each numbered model with useEM FALSE fits the unnumbered model given the
  # same variance components outside fast mode, by another solver, so the two
  # agree to rounding
  d = mvPop(nInd=40, nSnp=10)
  expect_warning(old <- RRBLUP2(d$pop, Vu=0.01, Ve=1, useEM=FALSE,
                                simParam=d$SP),
                 "RRBLUP2 is deprecated", class="deprecatedWarning")
  expect_equal(old, RRBLUP(d$pop, Vu=0.01, Ve=1, simParam=d$SP),
               tolerance=1e-8)
  expect_warning(old <- RRBLUP_D2(d$pop, Va=0.01, Vd=0.02, Ve=1,
                                  useEM=FALSE, simParam=d$SP),
                 "RRBLUP_D2 is deprecated", class="deprecatedWarning")
  expect_equal(old, RRBLUP_D(d$pop, Vu=c(0.01, 0.02), Ve=1,
                             simParam=d$SP), tolerance=1e-8)
  expect_warning(old <- RRBLUP_GCA2(d$pop, VuF=0.01, VuM=0.02, Ve=1,
                                    useEM=FALSE, simParam=d$SP),
                 "RRBLUP_GCA2 is deprecated", class="deprecatedWarning")
  expect_equal(old, RRBLUP_GCA(d$pop, Vu=c(0.01, 0.02), Ve=1,
                               simParam=d$SP), tolerance=1e-8)
  expect_warning(old <- RRBLUP_SCA2(d$pop, VuF=0.01, VuM=0.02, VuD=0.03,
                                    Ve=1, useEM=FALSE, simParam=d$SP),
                 "RRBLUP_SCA2 is deprecated", class="deprecatedWarning")
  expect_equal(old, RRBLUP_SCA(d$pop, Vu=c(0.01, 0.02, 0.03), Ve=1,
                               simParam=d$SP), tolerance=1e-8)
})

test_that("FIXED variance components from an earlier fit outside fast mode", {
  d = mvPop(nInd=40, nSnp=10)
  fit = RRBLUP(d$pop, simParam=d$SP)
  expect_identical(RRBLUP(d$pop, varComp=fit, simParam=d$SP),
                   RRBLUP(d$pop, Vu=fit@Vu, Ve=fit@Ve, simParam=d$SP))
  # Refitting with the REML estimates gives the REML fit back
  expect_equal(addEff(RRBLUP(d$pop, varComp=fit, simParam=d$SP)),
               addEff(fit), tolerance=1e-6)
  for(model in names(mkModels)){
    nTerm = if(model=="RRBLUP_SCA") 3L else 2L
    Vu = seq_len(nTerm)/100
    fit = mkModels[[model]](d$pop, Vu=Vu, Ve=1, simParam=d$SP)
    expect_equal(c(fit@Vu), Vu, info=model)
    # A fit's own Vu slot is passed back as it is
    expect_identical(mkModels[[model]](d$pop, varComp=fit, simParam=d$SP),
                     fit, info=model)
    expect_identical(mkModels[[model]](d$pop, Vu=fit@Vu, Ve=fit@Ve,
                                       simParam=d$SP),
                     fit, info=model)
    expect_error(mkModels[[model]](d$pop, Vu=Vu[-1], Ve=1, simParam=d$SP),
                 "one value per random term", info=model)
  }
})

test_that("FIXED variance components outside fast mode are for one trait", {
  # The direct solver for fixed variance components fits one trait, and fast
  # mode fits several
  d = mvPop(nInd=40, nSnp=10)
  fit = RRBLUP(d$pop, traits=1:2, simParam=d$SP)
  expect_error(RRBLUP(d$pop, traits=1:2, varComp=fit, simParam=d$SP),
               "only when fastMode is TRUE")
  expect_error(RRBLUP(d$pop, traits=1:2, Vu=fit@Vu, Ve=fit@Ve,
                      simParam=d$SP),
               "only when fastMode is TRUE")
  expect_error(RRBLUP_D(d$pop, traits=1:2, Vu=list(diag(2), diag(2)),
                        Ve=diag(2), simParam=d$SP),
               "only when fastMode is TRUE")
  expect_error(RRBLUPMemUse(nInd=100, nMarker=50, model="RRBLUP_D",
                            nTraits=2, useEM=FALSE),
               "only when fastMode is TRUE")
})

test_that("RRBLUPMemUse costs a model given its variance components as its REML fit", {
  # Held fixed, the variance components are evaluated once by the solver of
  # the multikernel REML fit, so its estimate is the bound used
  for(size in list(c(1000, 500), c(200, 5000))){
    for(model in c("RRBLUP_D","RRBLUP_GCA","RRBLUP_SCA")){
      expect_equal(RRBLUPMemUse(nInd=size[1], nMarker=size[2], model=model,
                                useEM=FALSE),
                   RRBLUPMemUse(nInd=size[1], nMarker=size[2], model=model),
                   info=paste(model, size[1]))
    }
    # RRBLUP uses that solver with a single term, which costs less than two
    v = RRBLUPMemUse(nInd=size[1], nMarker=size[2], model="RRBLUP",
                     useEM=FALSE)
    expect_gt(v, 0)
    expect_lt(v, RRBLUPMemUse(nInd=size[1], nMarker=size[2],
                              model="RRBLUP_D"))
  }
  # Fast mode does not take useEM
  expect_equal(RRBLUPMemUse(nInd=1000, nMarker=500, model="RRBLUP_D",
                            fastMode=TRUE, useEM=FALSE),
               RRBLUPMemUse(nInd=1000, nMarker=500, model="RRBLUP_D",
                            fastMode=TRUE))
})

test_that("DIPLOID the SCA models refuse other ploidy levels", {
  # The model splits each genome into one maternal and one paternal
  # haplotype, the design of a diploid single cross
  d4 = mvPop(nInd=40, nSnp=10, ploidy=4L)
  expect_error(RRBLUP_SCA(d4$pop, simParam=d4$SP), "only diploids")
  expect_error(RRBLUP_SCA(d4$pop, Vu=c(0.01, 0.01, 0.01), Ve=1,
                          simParam=d4$SP), "only diploids")
  expect_error(quietDeprecated(RRBLUP_SCA2(d4$pop, simParam=d4$SP)),
               "only diploids")
})

test_that("FIXED a variance of zero outside fast mode leaves its term out", {
  # REML often estimates a dominance or SCA variance at exactly zero, and a
  # fit like that has to be usable as varComp. A term with no variance has
  # no effects, and the other terms are fitted as if it were not there.
  d = mvPop(nInd=40, nSnp=10)
  dom = RRBLUP_D(d$pop, Vu=c(0.01, 0), Ve=1, simParam=d$SP)
  # A marker's dominance effect is then the directional dominance alone, the
  # same at every segregating marker
  domEff = dom@gv[[1]]@domEff
  expect_lte(length(unique(signif(domEff[domEff!=0], 10))), 1L)
  expect_true(all(is.finite(dom@gv[[1]]@addEff)))
  expect_identical(RRBLUP_D(d$pop, varComp=dom, simParam=d$SP), dom)
  # Without the maternal term only the paternal effects are fitted, which is
  # RRBLUP on the paternal haplotypes
  gca = RRBLUP_GCA(d$pop, Vu=c(0, 0.01), Ve=1, simParam=d$SP)
  expect_true(all(gca@gv[[1]]@addEff==0))
  M2 = 2*pullSnpHaplo(d$pop, haplo=2, simParam=d$SP) - 1
  y = d$pop@pheno[, 1L, drop=FALSE]
  # mvPop has three fixed effect levels
  X = model.matrix(~factor(d$pop@fixEff)-1)
  ref = solveRRBLUP_EM(y, X, M2, Vu=0.01, Ve=1, tol=1e-6, maxIter=10L,
                       useEM=FALSE)
  expect_equal(c(gca@gv[[1]]@addEffMale), c(ref$u), tolerance=1e-8)
  # Every term at zero leaves only the fixed effects
  sca = RRBLUP_SCA(d$pop, Vu=c(0, 0, 0), Ve=1, simParam=d$SP)
  expect_true(all(c(sca@gv[[1]]@addEff, sca@gv[[1]]@addEffMale)==0))
  domEff = sca@gv[[1]]@domEff
  expect_lte(length(unique(signif(domEff[domEff!=0], 10))), 1L)
  expect_true(all(RRBLUP(d$pop, Vu=0, Ve=1, simParam=d$SP)@gv[[1]]@addEff==0))
})
