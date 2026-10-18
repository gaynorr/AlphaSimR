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

test_that("RRBLUP takes variance components only in fast mode", {
  # The direct solver always estimates them, so supplying them without
  # fastMode would otherwise be ignored without notice
  d = gsPop(nInd=30, nQtl=5, nSnp=5)
  expect_error(RRBLUP(d$pop, Vu=0.01, simParam=d$SP), "fastMode")
  expect_error(RRBLUP(d$pop, Ve=1, simParam=d$SP), "fastMode")
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
  oneB = RRBLUP2(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1], useEM=FALSE,
                 nThreads=1L, simParam=d$SP)
  manyB = RRBLUP2(d$pop, Vu=ref@Vu[1,1], Ve=ref@Ve[1,1], useEM=FALSE,
                  nThreads=nThreads, simParam=d$SP)
  expect_equal(addEff(oneB), addEff(manyB), tolerance=1e-8)
})

# The models read their markers straight from the population's packed
# genotypes. Building the same marker matrices in R from the pull functions
# and handing them to the exported solvers is an independent route to the
# same fit, and with a single fixed effect level the design matrix is a
# column of ones on both routes. The codings below are exact in double
# precision, so the two routes have to agree to the last bit.

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
  ans2 = RRBLUP2(d$pop, Vu=0.1, Ve=1, useEM=FALSE, simParam=d$SP)
  expect_identical(addEff(ans2), c(ref2$u))
})

test_that("RRBLUP_D2 reads the same markers that pullSnpGeno returns", {
  d = gsPop(nInd=60, nQtl=8, nSnp=12, nChr=3)
  y = d$pop@pheno[, 1L, drop=FALSE]
  G = pullSnpGeno(d$pop, simParam=d$SP)
  Ma = G - 1
  Md = domCode(G, 2)
  X = cbind(1, rowSums(Md))
  ref = solveRRBLUP_EM2(y, X, Ma, Md, Vu1=0.1, Vu2=0.05, Ve=1, tol=1e-6,
                        maxIter=10L, useEM=FALSE)
  ans = RRBLUP_D2(d$pop, Va=0.1, Vd=0.05, Ve=1, useEM=FALSE, simParam=d$SP)
  expect_identical(addEff(ans), c(ref$u[, 1L]))
})

test_that("RRBLUP_GCA2 and RRBLUP_SCA2 read the parental haplotypes", {
  d = gsPop(nInd=60, nQtl=8, nSnp=12, nChr=3)
  y = d$pop@pheno[, 1L, drop=FALSE]
  X1 = matrix(1, nrow(y), 1L)
  M1 = 2*pullSnpHaplo(d$pop, haplo=1, simParam=d$SP) - 1
  M2 = 2*pullSnpHaplo(d$pop, haplo=2, simParam=d$SP) - 1
  ref = solveRRBLUP_EM2(y, X1, M1, M2, Vu1=0.1, Vu2=0.2, Ve=1, tol=1e-6,
                        maxIter=10L, useEM=FALSE)
  gca = RRBLUP_GCA2(d$pop, VuF=0.1, VuM=0.2, Ve=1, useEM=FALSE,
                    simParam=d$SP)
  expect_identical(c(gca@gv[[1L]]@addEff), c(ref$u[, 1L]))
  expect_identical(c(gca@gv[[1L]]@addEffMale), c(ref$u[, 2L]))

  M3 = domCode(pullSnpGeno(d$pop, simParam=d$SP), 2)
  ref3 = solveRRBLUP_EM3(y, cbind(1, rowSums(M3)), M1, M2, M3, Vu1=0.1,
                         Vu2=0.2, Vu3=0.05, Ve=1, tol=1e-6, maxIter=10L,
                         useEM=FALSE)
  sca = RRBLUP_SCA2(d$pop, VuF=0.1, VuM=0.2, VuD=0.05, Ve=1, useEM=FALSE,
                    simParam=d$SP)
  expect_identical(c(sca@gv[[1L]]@addEff), c(ref3$u[, 1L]))
  expect_identical(c(sca@gv[[1L]]@addEffMale), c(ref3$u[, 2L]))
})

test_that("RRBLUP_GCA2 splits a tetraploid's haplotypes in two", {
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
  gca = RRBLUP_GCA2(pop, VuF=0.1, VuM=0.2, Ve=1, useEM=FALSE, simParam=SP)
  expect_identical(c(gca@gv[[1L]]@addEff), c(ref$u[, 1L]))
  expect_identical(c(gca@gv[[1L]]@addEffMale), c(ref$u[, 2L]))
})

test_that("setEBV reads the parental haplotypes of GCA and SCA solutions", {
  skip_on_cran()
  d = gsPop(nInd=60, nQtl=8, nSnp=12, nChr=3)
  H1 = pullSnpHaplo(d$pop, haplo=1, simParam=d$SP)
  H2 = pullSnpHaplo(d$pop, haplo=2, simParam=d$SP)
  het = domCode(H1 + H2, 2)
  gcaSol = RRBLUP_GCA2(d$pop, VuF=0.1, VuM=0.2, Ve=1, useEM=FALSE,
                       simParam=d$SP)
  gca = gcaSol@gv[[1L]]
  expected = gca@intercept + (2*H1 - 1)%*%gca@addEff +
    (2*H2 - 1)%*%gca@addEffMale
  expect_equal(c(ebv(setEBV(d$pop, gcaSol, simParam=d$SP))), c(expected),
               tolerance=1e-10)
  scaSol = RRBLUP_SCA2(d$pop, VuF=0.1, VuM=0.2, VuD=0.05, Ve=1,
                       useEM=FALSE, simParam=d$SP)
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
    expect_error(fit(d$pop, Vu=rep(0.01, nTerm), simParam=d$SP),
                 "fastMode", info=model)
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
    ans = mkModels[[model]](d$pop, traits=1:2, fastMode=TRUE, Vu=Vu, Ve=Ve,
                            maxIter=10000L, simParam=d$SP)
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
  ans = RRBLUP_GCA(dh, fastMode=TRUE, Vu=c(0.02, 0.02), Ve=1e-16,
                   maxIter=10000L, simParam=d$SP)
  female = ans@female[[1]]@addEff
  male = ans@male[[1]]@addEff
  expect_true(all(is.finite(female)))
  expect_gt(max(abs(female)), 0)
  expect_equal(female, male, tolerance=1e-6)
  ref = RRBLUP(dh, fastMode=TRUE, Vu=0.04, Ve=1e-16, maxIter=10000L,
               simParam=d$SP)
  expect_equal(2*female, ref@bv[[1]]@addEff, tolerance=1e-6)
})
