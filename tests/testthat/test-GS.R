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
  # And fastRRBLUP's variance components from a record cross product
  expect_gt(RRBLUPMemUse(nInd=10000, nMarker=10000, model="fastRRBLUP",
                         nTraits=3, nSubsample=2000),
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

# Multiple trait fastRRBLUP. With the variance components supplied it is the
# BLUP of the multiple trait model, which Henderson's equations give directly
# for a population small enough to write them out in R. With every record
# used for the variance components it is the same REML fit as multiple trait
# RRBLUP, so the two have to agree.

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
# with the markers centered as fastRRBLUP centers them
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

test_that("multiple trait fastRRBLUP is the BLUP of the multiple trait model", {
  d = mvPop(nInd=80, nSnp=15)
  ans = fastRRBLUP(d$pop, traits=1:3, Vu=mvVu, Ve=mvVe, simParam=d$SP)
  ref = denseMvBlup(d$pop, d$SP, mvVu, mvVe)
  # It iterates rather than solving, so it is held to a looser tolerance
  expect_equal(unname(addEffAll(ans)), ref, tolerance=1e-6)
  expect_equal(length(ans@bv), 3L)
  expect_equal(dimnames(ans@Vu), list(d$SP$traitNames, d$SP$traitNames))
})

test_that("multiple trait fastRRBLUP reads a tetraploid's markers", {
  # A ploidy other than two goes through the general loop of the products
  d = mvPop(nInd=60, nSnp=12, ploidy=4L)
  ans = fastRRBLUP(d$pop, traits=1:3, Vu=mvVu, Ve=mvVe, simParam=d$SP)
  ref = denseMvBlup(d$pop, d$SP, mvVu, mvVe)
  expect_equal(unname(addEffAll(ans)), ref, tolerance=1e-6)
})

test_that("multiple trait fastRRBLUP agrees with RRBLUP on every record", {
  skip_on_cran()
  # More records than markers decomposes the markers, and fewer decomposes
  # the record cross product built from the genotypes
  for(size in list(c(nInd=80, nSnp=15), c(nInd=40, nSnp=30))){
    d = mvPop(nInd=size[["nInd"]], nSnp=size[["nSnp"]])
    ans = fastRRBLUP(d$pop, traits=1:3, nSubsample=0L, simParam=d$SP)
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

test_that("fastRRBLUP predicts the same genetic values as RRBLUP", {
  # The intercept of the genetic values goes with the uncentered markers.
  # Taking it from the centered ones shifts every prediction by Mmean*u,
  # which is not zero in general, so this checks a trait where it is not.
  d = mvPop(nInd=80, nSnp=15)
  ans = fastRRBLUP(d$pop, traits=1, nSubsample=0L, simParam=d$SP)
  ref = RRBLUP(d$pop, traits=1, simParam=d$SP)
  expect_gt(abs(ref@bv[[1]]@intercept), 1e-8)
  expect_equal(ans@gv[[1]]@intercept, ref@gv[[1]]@intercept, tolerance=1e-6)
  expect_equal(c(ebv(setEBV(d$pop, ans, simParam=d$SP))),
               c(ebv(setEBV(d$pop, ref, simParam=d$SP))), tolerance=1e-6)
})

test_that("multiple trait fastRRBLUP leaves its inputs unchanged", {
  d = mvPop(nInd=40, nSnp=10)
  Vu = mvVu + 0
  Ve = mvVe + 0
  pheno0 = d$pop@pheno + 0
  fastRRBLUP(d$pop, traits=1:3, Vu=Vu, Ve=Ve, simParam=d$SP)
  expect_identical(Vu, mvVu)
  expect_identical(Ve, mvVe)
  expect_identical(d$pop@pheno, pheno0)
})

test_that("multiple trait fastRRBLUP rejects unusable variance components", {
  d = mvPop(nInd=40, nSnp=10)
  expect_error(fastRRBLUP(d$pop, traits=1:3, Vu=diag(2), Ve=diag(3),
                          simParam=d$SP), "one row per trait")
  expect_error(fastRRBLUP(d$pop, traits=1:3, Vu=mvVu, Ve=diag(c(1,1,0)),
                          simParam=d$SP), "positive definite")
  expect_error(fastRRBLUP(d$pop, traits=1:3, Vu=-mvVu, Ve=mvVe,
                          simParam=d$SP), "positive semidefinite")
  expect_error(fastRRBLUP(d$pop, traits=1:3, Vu=0*mvVu, Ve=mvVe,
                          simParam=d$SP), "not positive")
})

test_that("multiple trait fastRRBLUP rejects an indefinite Vu", {
  d = mvPop(nInd=40, nSnp=10)
  expect_error(fastRRBLUP(d$pop, traits=1:3, Vu=diag(c(0.01,0.01,-0.005)),
                          Ve=mvVe, simParam=d$SP), "positive semidefinite")
})

test_that("multiple trait fastRRBLUP checks Vu on its own scale", {
  # An error variance close to zero for the second trait would make the
  # first trait's negative variance look like rounding after whitening
  d = mvPop(nInd=40, nSnp=10)
  expect_error(fastRRBLUP(d$pop, traits=1:2, Vu=diag(c(-0.01, 0.01)),
                          Ve=diag(c(1, 1e-14)), simParam=d$SP),
               "positive semidefinite")
})

test_that("multiple trait fastRRBLUP warns when maxIter is reached", {
  d = mvPop(nInd=80, nSnp=15)
  expect_warning(fastRRBLUP(d$pop, traits=1:3, Vu=mvVu, Ve=mvVe, maxIter=1,
                            simParam=d$SP), "without converging")
})

test_that("multiple trait fastRRBLUP with diagonal covariance matrices fits each trait on its own", {
  # With no covariance between the traits, the multiple trait BLUP of each
  # trait is its single trait BLUP. An error variance close to zero for the
  # first trait makes its variance ratio far larger than the second's, which
  # must not stop the second trait from being fitted or from converging.
  d = mvPop(nInd=80, nSnp=15)
  Vu = diag(c(0.02, 0.03))
  Ve = diag(c(1e-14, 1))
  ans = fastRRBLUP(d$pop, traits=1:2, Vu=Vu, Ve=Ve, maxIter=10000,
                   simParam=d$SP)
  one = fastRRBLUP(d$pop, traits=2, Vu=Vu[2,2], Ve=Ve[2,2], maxIter=10000,
                   simParam=d$SP)
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
  ans = fastRRBLUP(d$pop, traits=pair, Vu=mvVu[1:2, 1:2], Ve=mvVe[1:2, 1:2],
                   simParam=d$SP)
  expect_equal(dimnames(ans@Vu), list(c("Response1","Response2"),
                                      c("Response1","Response2")))
  expect_equal(ans@gv[[2]]@name, "est_GV_Response2")
  ans = RRBLUP_D(d$pop, traits=pair, simParam=d$SP)
  expect_equal(rownames(ans@Vu)[1:2], c("a:Response1", "a:Response2"))
})

test_that("models accept a function of the traits that returns a vector", {
  d = mvPop(nInd=60, nSnp=12)
  index = function(x) rowSums(x)
  ans = fastRRBLUP(d$pop, traits=index, Vu=0.01, Ve=1, simParam=d$SP)
  expect_equal(length(ans@gv), 1L)
  ans = fastRRBLUP(d$pop, traits=index, simParam=d$SP)
  expect_true(all(is.finite(addEff(ans))))
  for(fit in list(RRBLUP_D, RRBLUP_GCA, RRBLUP_SCA)){
    ans = fit(d$pop, traits=index, simParam=d$SP)
    expect_equal(length(ans@gv), 1L)
    expect_equal(dim(ans@Ve), c(1L, 1L))
  }
})

test_that("multiple trait fastRRBLUP needs every trait recorded", {
  expect_error(AlphaSimR:::checkCompleteRecords(matrix(c(1,NA,3,4), 2)),
               "every trait recorded")
  expect_null(AlphaSimR:::checkCompleteRecords(matrix(1:4, 2)))
})

test_that("multiple trait fastRRBLUP is not affected by the number of threads", {
  skip_on_cran()
  skip_if(getNumThreads() < 2L)
  d = mvPop(nInd=80, nSnp=15)
  one = fastRRBLUP(d$pop, traits=1:3, nThreads=1L, simParam=d$SP)
  two = fastRRBLUP(d$pop, traits=1:3, nThreads=2L, simParam=d$SP)
  expect_identical(addEffAll(one), addEffAll(two))
  expect_identical(one@Vu, two@Vu)
})
