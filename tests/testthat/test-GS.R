context("GS")

# The genomic selection functions had no tests at all. These do not check the
# models against an outside implementation, because there isn't one here.
# They use the redundancy inside the module instead: the same model can be
# reached by more than one route, and the routes have to agree.
#
# The most important of those is RRBLUP2. It solves the mixed model in one of
# two ways, working with a square matrix of the records when they are fewer
# than the markers and of the markers otherwise, and the pair only means
# anything if both give the same answer as RRBLUP does.

# A population with more records than markers, so that RRBLUP2 works with the
# markers, and one with fewer, so that it works with the records
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

test_that("RRBLUP2 agrees with RRBLUP when it works with the markers", {
  # 40 records against 10 markers, so Henderson's equations are cheaper
  d = gsPop(nInd=40, nQtl=5, nSnp=5)
  expect_lt(d$nMarker, nInd(d$pop))

  ref = RRBLUP(d$pop, simParam=d$SP)
  ans = RRBLUP2(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1], useEM=FALSE,
                simParam=d$SP)

  expect_equal(addEff(ans), addEff(ref), tolerance=1e-6)
  expect_equal(ans@gv[[1]]@intercept, ref@gv[[1]]@intercept, tolerance=1e-6)
})

test_that("RRBLUP2 agrees with RRBLUP when it works with the records", {
  # 15 records against 40 markers, so the record side is cheaper. This is the
  # branch that never forms a markers by records matrix, and the one a
  # mistake in the order of that product would show up in.
  d = gsPop(nInd=15, nQtl=5, nSnp=20)
  expect_gt(d$nMarker, nInd(d$pop))

  ref = RRBLUP(d$pop, simParam=d$SP)
  ans = RRBLUP2(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1], useEM=FALSE,
                simParam=d$SP)

  expect_equal(addEff(ans), addEff(ref), tolerance=1e-6)
  expect_equal(ans@gv[[1]]@intercept, ref@gv[[1]]@intercept, tolerance=1e-6)
})

test_that("both RRBLUP2 methods handle more than one fixed effect level", {
  # The record side splits the solution of V into the fixed effect columns
  # and the response, so it has to keep them apart when there is more than
  # one fixed effect. Both shapes of the problem are checked.
  for(size in list(c(nInd=40, nQtl=5, nSnp=5),
                   c(nInd=15, nQtl=5, nSnp=20))){
    d = gsPop(nInd=size[["nInd"]], nQtl=size[["nQtl"]], nSnp=size[["nSnp"]],
              nFixEff=3L)
    expect_equal(length(unique(d$pop@fixEff)), 3L)

    ref = RRBLUP(d$pop, simParam=d$SP)
    ans = RRBLUP2(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1], useEM=FALSE,
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

  ans = RRBLUP2(d$pop, useEM=TRUE, maxIter=200, simParam=d$SP)

  expect_true(is.finite(ans@Vu[1,1]))
  expect_true(is.finite(ans@Ve[1,1]))
  expect_gt(ans@Vu[1,1], 0)
  expect_gt(ans@Ve[1,1], 0)
  expect_true(all(is.finite(addEff(ans))))

  # Predictions from an estimated model should still track the truth
  pred = setEBV(d$pop, ans, simParam=d$SP)
  expect_gt(cor(c(ebv(pred)), c(gv(d$pop))), 0.3)
})

test_that("RRBLUP2 gives the same answer whatever its starting values", {
  # Supplying the variance components turns estimation off, so two calls with
  # the same components must agree exactly and a third with different
  # components must not
  d = gsPop(nInd=30, nQtl=5, nSnp=5)
  ref = RRBLUP(d$pop, simParam=d$SP)

  a = RRBLUP2(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1], useEM=FALSE,
              simParam=d$SP)
  b = RRBLUP2(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1], useEM=FALSE,
              simParam=d$SP)
  expect_identical(addEff(a), addEff(b))

  c2 = RRBLUP2(d$pop, Vu=ref@Vu[1,1]*10, Ve=ref@Ve[1,1], useEM=FALSE,
               simParam=d$SP)
  expect_false(isTRUE(all.equal(addEff(a), addEff(c2))))

  # A larger marker variance shrinks the effects less
  expect_gt(sum(addEff(c2)^2), sum(addEff(a)^2))
})

test_that("fastRRBLUP agrees with RRBLUP given the same variance components", {
  # Iterates to convergence rather than solving, so it is the slowest fit here
  skip_on_cran()
  d = gsPop(nInd=40, nQtl=5, nSnp=5)
  ref = RRBLUP(d$pop, simParam=d$SP)
  ans = fastRRBLUP(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1], maxIter=10000,
                   simParam=d$SP)

  # fastRRBLUP iterates rather than solving directly, so it is held to a
  # looser tolerance than the direct solvers are held to each other
  expect_equal(addEff(ans), addEff(ref), tolerance=1e-4)
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
  # dosages are centred before they are weighted rather than used raw.
  M = pullSnpGeno(d$pop, simParam=d$SP)
  a = addEff(ans)
  p = d$pop@ploidy
  centred = (M - p/2) * (2/p)
  byHand = c(centred %*% a) + ans@gv[[1]]@intercept
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

test_that("the numbered dominance and hybrid models run with fixed components", {
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

  dom = RRBLUP_D2(pop, Va=vg/20, Vd=vg/20, Ve=ve, useEM=FALSE, simParam=SP)
  expect_true(all(is.finite(addEff(dom))))

  gca = RRBLUP_GCA2(pop, VuF=vg/20, VuM=vg/20, Ve=ve, useEM=FALSE,
                    simParam=SP)
  expect_true(all(is.finite(addEff(gca))))

  sca = RRBLUP_SCA2(pop, VuF=vg/20, VuM=vg/20, VuD=vg/20, Ve=ve,
                    useEM=FALSE, simParam=SP)
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

  expect_error(RRBLUPMemUse(nInd=100, nMarker=50, model="notAModel"),
               "not recognized")

  # Every estimate grows with the size of the problem
  for(m in models){
    small = RRBLUPMemUse(nInd=500, nMarker=500, model=m)
    large = RRBLUPMemUse(nInd=2000, nMarker=2000, model=m)
    expect_gt(large, small)
  }
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
  oneB = RRBLUP2(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1], useEM=FALSE,
                 nThreads=1L, simParam=d$SP)
  manyB = RRBLUP2(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1], useEM=FALSE,
                  nThreads=nThreads, simParam=d$SP)
  expect_equal(addEff(oneB), addEff(manyB), tolerance=1e-8)
})
