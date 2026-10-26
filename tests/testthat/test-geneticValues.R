context("geneticValues")

# Genetic values worked out in R from the QTL dosages, using the coding the
# trait classes document. These are the values the C++ code must reproduce
# whichever method it picks for a trait.
refGv = function(pop, simParam, trait) {
  tr = simParam$traits[[trait]]
  X = pullQtlGeno(pop, trait = trait, simParam = simParam)
  P = pop@ploidy
  xa = (X - P / 2) * (2 / P)
  xd = X * (P - X) * (2 / P)^2
  gv = tr@intercept + drop(xa %*% tr@addEff)
  if (.hasSlot(tr, "domEff")) {
    gv = gv + drop(xd %*% tr@domEff)
  }
  if (.hasSlot(tr, "epiEff")) {
    E = tr@epiEff
    gv = gv + drop((xa[, E[, 1], drop = FALSE] * xa[, E[, 2], drop = FALSE]) %*% E[, 3])
  }
  gxe = NULL
  if (.hasSlot(tr, "gxeEff")) {
    gxe = unname(tr@gxeInt + drop(xa %*% tr@gxeEff))
  }
  list(gv = unname(gv), gxe = gxe)
}

expectMatchesRef = function(pop, simParam) {
  for (i in seq_len(simParam$nTraits)) {
    ref = refGv(pop, simParam, i)
    expect_equal(unname(pop@gv[, i]), ref$gv)
    if (!is.null(ref$gxe)) {
      expect_equal(as.vector(pop@gxe[[i]]), ref$gxe)
    }
  }
}

# ---- 1. Agreement with the dosage coding ----

test_that("MATCH dense additive, dominance and GxE traits", {
  # Every site is a QTL, so most lookup tables cover several loci
  set.seed(101)
  founderPop = quickHaplo(nInd = 20, nChr = 3, segSites = 64)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(64)
  SP$addTraitAD(64, meanDD = 0.3, varDD = 0.2)
  SP$addTraitADG(64, meanDD = 0.2, varGxE = 0.5)
  pop = randCross(newPop(founderPop, simParam = SP), 50, simParam = SP)
  expectMatchesRef(pop, SP)
})

test_that("MATCH correlated traits sharing dense and sparse loci", {
  # Five traits added together share their loci and are worked out in one
  # pass. With every site a QTL the pass uses lookup tables, and with one site
  # in ten it uses dosages, so both methods are checked for a group.
  set.seed(102)
  founderPop = quickHaplo(nInd = 20, nChr = 2, segSites = 200)
  corA = diag(5) * 0.5 + 0.5
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(200, mean = rep(0, 5), var = rep(1, 5), corA = corA)
  SP$addTraitADG(20, mean = rep(0, 5), var = rep(1, 5), meanDD = rep(0.2, 5),
                 varGxE = rep(0.3, 5), corA = corA)
  pop = randCross(newPop(founderPop, simParam = SP), 50, simParam = SP)
  expectMatchesRef(pop, SP)
})

test_that("MATCH epistatic traits", {
  set.seed(103)
  founderPop = quickHaplo(nInd = 20, nChr = 3, segSites = 40)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitAE(20, relAA = 0.5)
  SP$addTraitADEG(20, meanDD = 0.2, relAA = 0.5, varGxE = 0.3)
  pop = randCross(newPop(founderPop, simParam = SP), 50, simParam = SP)
  expectMatchesRef(pop, SP)
})

test_that("MATCH epistatic traits on sparse loci", {
  # Few loci share a byte here, so pairs read the genotypes directly rather
  # than through a buffer of dosages, which the test above uses
  set.seed(106)
  founderPop = quickHaplo(nInd = 20, nChr = 3, segSites = 400)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitADEG(20, meanDD = 0.2, relAA = 0.5, varGxE = 0.3)
  pop = randCross(newPop(founderPop, simParam = SP), 50, simParam = SP)
  expectMatchesRef(pop, SP)
})

test_that("MATCH polyploid traits", {
  # Additive traits use lookup tables at any ploidy, while dominance above
  # diploid can only use dosages
  set.seed(104)
  founderPop = quickHaplo(nInd = 20, nChr = 2, segSites = 40, ploidy = 4)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(40)
  SP$addTraitAD(40, meanDD = 0.3, varDD = 0.2)
  SP$addTraitAE(20, relAA = 0.5)
  pop = randCross(newPop(founderPop, simParam = SP), 30, simParam = SP)
  expectMatchesRef(pop, SP)
})

test_that("MATCH chromosomes that do not fill their last byte", {
  set.seed(105)
  founderPop = quickHaplo(nInd = 5, nChr = 3, segSites = 13)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitAD(5, meanDD = 0.3)
  SP$addTraitA(1)
  pop = randCross(newPop(founderPop, simParam = SP), 7, simParam = SP)
  expectMatchesRef(pop, SP)
  expectMatchesRef(pop[1], SP)
})

# ---- 2. Consistency between paths ----

test_that("AGREE resetPop reproduces the values newPop gives", {
  # Correlated traits are worked out together, and resetPop has to make the
  # same choice of method for them to give identical values. Both groups
  # here are sparse enough that a group uses dosages while any one of its
  # traits on its own would use lookup tables, so evaluating trait by trait
  # would give different values in the last few bits.
  set.seed(201)
  founderPop = quickHaplo(nInd = 20, nChr = 2, segSites = 200)
  corA = matrix(c(1, 0.5, 0.3, 0.5, 1, 0.2, 0.3, 0.2, 1), 3)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(20, mean = c(0, 0, 0), var = c(1, 1, 1), corA = corA)
  SP$addTraitADG(100, mean = c(0, 1, 2), var = c(1, 2, 1),
                 meanDD = c(0.1, 0.3, 0), varGxE = c(0.3, 0.5, 0.2),
                 corA = corA)
  SP$addTraitA(200)
  SP$addTraitAE(20, relAA = 0.5)
  pop = randCross(newPop(founderPop, simParam = SP), 50, simParam = SP)
  reset = resetPop(pop, simParam = SP)
  expect_identical(gv(reset), gv(pop))
  expect_identical(reset@gxe, pop@gxe)
})

test_that("THREADS genetic values do not depend on the thread count", {
  skip_if(getNumThreads() < 2)
  values = lapply(c(1L, 2L), function(nThreads) {
    set.seed(202)
    founderPop = quickHaplo(nInd = 20, nChr = 2, segSites = 200)
    SP = SimParam$new(founderPop)
    SP$nThreads = nThreads
    SP$addTraitA(200, mean = c(0, 0), var = c(1, 1),
                 corA = matrix(c(1, 0.5, 0.5, 1), 2))
    SP$addTraitADG(20, meanDD = 0.2, varGxE = 0.3)
    SP$addTraitADE(20, meanDD = 0.2, relAA = 0.5)
    pop = randCross(newPop(founderPop, simParam = SP), 200, simParam = SP)
    list(gv = gv(pop), gxe = pop@gxe)
  })
  expect_identical(values[[1]], values[[2]])
})

# ---- 3. Genetic parameters ----

# Breeding values and dominance deviations worked out in R from the QTL
# dosages and the observed genotype frequencies at each locus, following the
# definitions calcGenParam documents
refBvDd = function(pop, simParam, trait) {
  tr = simParam$traits[[trait]]
  X = unname(pullQtlGeno(pop, trait = trait, simParam = simParam))
  P = pop@ploidy
  x = 0:P
  xa = (x - P / 2) * (2 / P)
  xd = x * (P - x) * (2 / P)^2
  d = if (.hasSlot(tr, "domEff")) tr@domEff else rep(0, ncol(X))
  bv = dd = numeric(nrow(X))
  alpha = numeric(ncol(X))
  for (i in seq_len(ncol(X))) {
    freq = tabulate(X[, i] + 1, nbins = P + 1) / nrow(X)
    genoMu = sum(freq * x)
    g = xa * tr@addEff[i] + xd * d[i]
    gvMu = sum(freq * g)
    alpha[i] = sum(freq * (g - gvMu) * (x - genoMu)) / sum(freq * (x - genoMu)^2)
    if (!is.finite(alpha[i])) alpha[i] = 0
    b = (x - genoMu) * alpha[i]
    bv = bv + b[X[, i] + 1]
    dd = dd + (g - b - gvMu)[X[, i] + 1]
  }
  list(bv = bv, dd = dd, alpha = alpha)
}

# Minimum norm least squares on columns scaled to unit variance, leaving out
# columns with no variance, as calcGenParamE does when the dosages do not
# identify separate slopes. Equal weights for the individuals give the
# observed joint genotype frequencies.
refMinNormFit = function(Z, y) {
  coef = numeric(ncol(Z))
  ss = colMeans(Z^2)
  keep = which(ss > 1e-14)
  if (length(keep) == 0) return(coef)
  Zs = sweep(Z[, keep, drop = FALSE], 2, sqrt(ss[keep]), "/")
  e = eigen(crossprod(Zs) / nrow(Zs), symmetric = TRUE)
  use = e$values > 1e-10 * max(e$values)
  b = e$vectors[, use, drop = FALSE] %*%
    (crossprod(e$vectors[, use, drop = FALSE], colMeans(Zs * y)) / e$values[use])
  coef[keep] = as.vector(b) / sqrt(ss[keep])
  coef
}

# The same for a trait with epistasis, following calcGenParamE: each locus's
# own effects are regressed on its own dosage, each pair's epistatic value
# is fitted to both dosages jointly for its part of the breeding values, and
# to the dosages and dominance codings for its part of the dominance
# deviations, and the epistatic deviation is what is left of the pair's
# genetic value
refBvDdAa = function(pop, simParam, trait) {
  tr = simParam$traits[[trait]]
  X = unname(pullQtlGeno(pop, trait = trait, simParam = simParam))
  P = pop@ploidy
  xa = function(x) (x - P / 2) * (2 / P)
  xd = function(x) x * (P - x) * (2 / P)^2
  d = if (.hasSlot(tr, "domEff")) tr@domEff else rep(0, ncol(X))
  E = tr@epiEff
  n = nrow(X)
  bv = dd = aa = numeric(n)
  alpha = numeric(ncol(X))
  for (i in seq_len(nrow(E))) {
    l = E[i, 1:2]
    Zx = sweep(X[, l, drop = FALSE], 2, colMeans(X[, l, drop = FALSE]))
    Zd = xd(X[, l, drop = FALSE])
    Zd = sweep(Zd, 2, colMeans(Zd))
    for (k in 1:2) {
      own = xa(X[, l[k]]) * tr@addEff[l[k]] + xd(X[, l[k]]) * d[l[k]]
      own = own - mean(own)
      al = sum(own * Zx[, k]) / sum(Zx[, k]^2)
      if (!is.finite(al)) al = 0
      alpha[l[k]] = al
      bv = bv + al * Zx[, k]
      dd = dd + own - al * Zx[, k]
    }
    # The own effects leave nothing in the epistatic deviation, so it is
    # what the two fits leave of the centered epistatic value
    y = E[i, 3] * xa(X[, l[1]]) * xa(X[, l[2]])
    y = y - mean(y)
    gamma = refMinNormFit(Zx, y)
    beta = refMinNormFit(cbind(Zx, Zd), y)
    alpha[l] = alpha[l] + gamma
    pairBv = as.vector(Zx %*% gamma)
    pairDd = as.vector(cbind(Zx, Zd) %*% beta) - pairBv
    bv = bv + pairBv
    dd = dd + pairDd
    aa = aa + y - pairBv - pairDd
  }
  list(bv = bv, dd = dd, aa = aa, alpha = alpha)
}

expectGenParamMatchesRef = function(pop, simParam) {
  gp = genParam(pop, simParam = simParam)
  for (i in seq_len(simParam$nTraits)) {
    if (.hasSlot(simParam$traits[[i]], "epiEff")) {
      ref = refBvDdAa(pop, simParam, i)
      expect_equal(unname(gp$aa[, i]), ref$aa)
    } else {
      ref = refBvDd(pop, simParam, i)
    }
    expect_equal(unname(gp$bv[, i]), ref$bv)
    expect_equal(unname(gp$dd[, i]), ref$dd)
    expect_equal(as.vector(gp$alpha[[i]]), ref$alpha)
  }
}

test_that("MATCH breeding values and dominance deviations", {
  # Dense and sparse traits, and a correlated group, so lookup tables and
  # dosages are both used, alone and for traits worked out together
  set.seed(301)
  founderPop = quickHaplo(nInd = 20, nChr = 2, segSites = 200)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitAD(200, meanDD = 0.3, varDD = 0.2)
  SP$addTraitA(20)
  SP$addTraitAD(20, mean = rep(0, 3), var = rep(1, 3), meanDD = c(0.1, 0.2, 0.3),
                corA = diag(3) * 0.5 + 0.5)
  pop = randCross(newPop(founderPop, simParam = SP), 100, simParam = SP)
  expectGenParamMatchesRef(pop, SP)
})

test_that("MATCH breeding values for polyploids and monomorphic loci", {
  set.seed(302)
  founderPop = quickHaplo(nInd = 20, nChr = 2, segSites = 40, ploidy = 4)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(40)
  SP$addTraitAD(40, meanDD = 0.3, varDD = 0.2)
  pop = randCross(newPop(founderPop, simParam = SP), 60, simParam = SP)
  expectGenParamMatchesRef(pop, SP)

  # With three founders, many loci are fixed among their progeny and have
  # an average effect of zero
  set.seed(303)
  founderPop = quickHaplo(nInd = 3, nChr = 2, segSites = 100)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitAD(100, meanDD = 0.3)
  pop = randCross(newPop(founderPop, simParam = SP), 20, simParam = SP)
  expectGenParamMatchesRef(pop, SP)
})

test_that("MATCH breeding values and deviations for epistatic traits", {
  # Dense enough for the dosage buffer, sparse enough for direct reads, a
  # tetraploid, and few founders so that some loci are fixed
  set.seed(306)
  founderPop = quickHaplo(nInd = 20, nChr = 2, segSites = 200)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitAE(200, relAA = 0.5)
  SP$addTraitADE(20, meanDD = 0.3, varDD = 0.2, relAA = 0.5)
  pop = randCross(newPop(founderPop, simParam = SP), 100, simParam = SP)
  expectGenParamMatchesRef(pop, SP)

  set.seed(307)
  founderPop = quickHaplo(nInd = 20, nChr = 2, segSites = 40, ploidy = 4)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitADE(20, meanDD = 0.3, varDD = 0.2, relAA = 0.5)
  pop = randCross(newPop(founderPop, simParam = SP), 60, simParam = SP)
  expectGenParamMatchesRef(pop, SP)

  set.seed(308)
  founderPop = quickHaplo(nInd = 3, nChr = 2, segSites = 100)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitADE(100, meanDD = 0.3, relAA = 0.5)
  pop = randCross(newPop(founderPop, simParam = SP), 20, simParam = SP)
  expectGenParamMatchesRef(pop, SP)
})

test_that("JOINT epistatic pairs in LD match least-squares projections", {
  # MaCS founders have LD between linked loci, and pairs are forced onto
  # neighboring loci so that every pair is in strong LD. For each pair, the
  # epistatic value's part of the breeding values is its least-squares fit
  # to both dosages, and its part of the dominance deviations is what adding
  # the dominance codings adds to that fit. lm finds both by QR, a different
  # route from calcGenParamE, and its fitted values are unique.
  set.seed(309)
  founderPop = runMacs(nInd = 300, nChr = 1, segSites = 40)
  for (type in c("AE", "ADE")) {
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    if (type == "AE") {
      SP$addTraitAE(40, relAA = 0.5)
    } else {
      SP$addTraitADE(40, meanDD = 0.3, varDD = 0.2, relAA = 0.5)
    }
    trait = SP$traits[[1]]
    trait@epiEff[, 1:2] = matrix(1:40, ncol = 2, byrow = TRUE)
    SP2 = SimParam$new(founderPop)
    SP2$nThreads = 1L
    SP2$manAddTrait(trait)
    pop = newPop(founderPop, simParam = SP2)
    gp = genParam(pop, simParam = SP2)
    X = pullQtlGeno(pop, simParam = SP2)
    xa = function(x) x - 1
    xd = function(x) x * (2 - x)
    d = if (type == "ADE") trait@domEff else rep(0, 40)
    bv = dd = aa = numeric(nInd(pop))
    for (i in seq_len(nrow(trait@epiEff))) {
      l = trait@epiEff[i, 1:2]
      x1 = X[, l[1]]
      x2 = X[, l[2]]
      for (k in 1:2) {
        x = X[, l[k]]
        own = trait@addEff[l[k]] * xa(x) + d[l[k]] * xd(x)
        fit = fitted(lm(own ~ x))
        bv = bv + fit - mean(own)
        dd = dd + own - fit
      }
      y = trait@epiEff[i, 3] * xa(x1) * xa(x2)
      fitA = fitted(lm(y ~ x1 + x2))
      fitAD = fitted(lm(y ~ x1 + x2 + xd(x1) + xd(x2)))
      bv = bv + fitA - mean(y)
      dd = dd + fitAD - fitA
      aa = aa + y - fitAD
    }
    expect_equal(unname(gp$bv[, 1]), unname(bv), tolerance = 1e-8)
    expect_equal(unname(gp$dd[, 1]), unname(dd), tolerance = 1e-8)
    expect_equal(unname(gp$aa[, 1]), unname(aa), tolerance = 1e-8)
    # Matching lm also means each pair's epistatic deviation is uncorrelated
    # with its own two dosages and dominance codings. LD with the loci of
    # other pairs is not removed, and appears in the covariances.
    expect_equal(unname(gp$gv[, 1]),
                 unname(gp$mu[1] + gp$bv[, 1] + gp$dd[, 1] + gp$aa[, 1]),
                 tolerance = 1e-10)
  }
})

test_that("LD dominance deviations without dominance effects are LD covariance", {
  # A trait without dominance effects has none under linkage equilibrium,
  # so all of its dominance variance comes from LD within pairs
  set.seed(310)
  founderPop = runMacs(nInd = 300, nChr = 1, segSites = 40)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitAE(40, relAA = 0.5)
  trait = SP$traits[[1]]
  trait@epiEff[, 1:2] = matrix(1:40, ncol = 2, byrow = TRUE)
  SP2 = SimParam$new(founderPop)
  SP2$nThreads = 1L
  SP2$manAddTrait(trait)
  pop = newPop(founderPop, simParam = SP2)
  gp = genParam(pop, simParam = SP2)
  expect_gt(gp$varD[1, 1], 0)
  expect_equal(unname(gp$genicVarD), 0)
  expect_equal(unname(gp$covD_HW), 0)
  expect_equal(unname(gp$covD_L), gp$varD[1, 1], tolerance = 1e-10)
  expect_equal(unname(gp$covD_LMat[1, 1]), gp$varD[1, 1], tolerance = 1e-10)
  expect_equal(unname(gp$nd[, 1]), unname(gp$dd[, 1] + gp$aa[, 1]), tolerance = 1e-10)
  expect_equal(unname(gp$gv_d[, 1]), rep(0, nInd(pop)))
})

test_that("MEAN splits into its Hardy-Weinberg and LD parts", {
  # mu - mu_LE is the change in the mean from LD, which only epistatic
  # effects feel: each pair's effect times the covariance of its additive
  # codings. mu_LE - mu_HW is the change from departures from HWE, which
  # only dominance effects feel. Selection and few founders give both.
  for (ploidy in c(2L, 4L)) {
    set.seed(314 + ploidy)
    founderPop = runMacs(nInd = 100, nChr = 1, segSites = 40, ploidy = ploidy)
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    SP$addTraitADE(40, meanDD = 0.3, varDD = 0.2, relAA = 0.5)
    pop = selectCross(newPop(founderPop, simParam = SP), nInd = 10,
                      nCrosses = 100, use = "gv", simParam = SP)
    gp = genParam(pop, simParam = SP)
    trait = SP$traits[[1]]
    X = pullQtlGeno(pop, simParam = SP)
    xa = (X - ploidy / 2) * (2 / ploidy)
    xd = X * (ploidy - X) * (2 / ploidy)^2
    p = colMeans(X) / ploidy
    # Expected dominance coding under HWE: binomial genotype frequencies
    xdHW = sapply(p, function(pp) {
      k = 0:ploidy
      sum(dbinom(k, ploidy, pp) * k * (ploidy - k) * (2 / ploidy)^2)
    })
    covXa = function(i, j) mean((xa[, i] - mean(xa[, i])) * (xa[, j] - mean(xa[, j])))
    epiLD = sum(apply(trait@epiEff, 1, function(e) e[3] * covXa(e[1], e[2])))
    domHW = sum(trait@domEff * (colMeans(xd) - xdHW))
    expect_equal(unname(gp$mu - gp$mu_LE), epiLD, tolerance = 1e-10)
    expect_equal(unname(gp$mu_LE - gp$mu_HW), domHW, tolerance = 1e-10)
  }
  # Without epistasis LD cannot change the mean
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitAD(40, meanDD = 0.3)
  pop = selectCross(newPop(founderPop, simParam = SP), nInd = 10,
                    nCrosses = 100, use = "gv", simParam = SP)
  gp = genParam(pop, simParam = SP)
  expect_equal(gp$mu_LE, gp$mu)
})

test_that("CONTINUITY a zero epistatic effect gives the trait without epistasis", {
  # With every epistatic effect zero, the joint fit has nothing to fit, so
  # the decomposition must equal that of the same additive and dominance
  # effects without epistasis, whatever the LD
  set.seed(311)
  founderPop = runMacs(nInd = 200, nChr = 1, segSites = 40)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitADE(40, meanDD = 0.3, varDD = 0.2, relAA = 0.5)
  pop = newPop(founderPop, simParam = SP)
  epi = SP$traits[[1]]
  epi@epiEff[, 3] = 0
  noEpi = as(epi, "TraitAD", strict = TRUE)
  withEpi = AlphaSimR:::calcGenParam(epi, pop, 1L)
  without = AlphaSimR:::calcGenParam(noEpi, pop, 1L)
  expect_equal(as.vector(withEpi$bv), as.vector(without$bv), tolerance = 1e-10)
  expect_equal(as.vector(withEpi$dd), as.vector(without$dd), tolerance = 1e-10)
  expect_equal(as.vector(withEpi$alpha), as.vector(without$alpha), tolerance = 1e-10)
  expect_equal(as.vector(withEpi$aa), rep(0, nInd(pop)), tolerance = 1e-10)
})

test_that("SYMMETRY perfect LD splits the slope the same way in either order", {
  # The second locus of each pair copies the first, or its complement, so
  # only the pair's total breeding value is identified. The split between
  # the loci must not depend on which is listed first.
  set.seed(312)
  h1 = rbinom(400, 1, 0.6)
  h3 = rbinom(400, 1, 0.3)
  haplo = cbind(h1, h1, h3, 1 - h3)
  colnames(haplo) = paste0("m", 1:4)
  genMap = data.frame(marker = colnames(haplo), chr = 1, pos = c(0.1, 0.2, 0.3, 0.4))
  founderPop = importHaplo(haplo, genMap)
  fit = function(pairs) {
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    SP$importTrait(markerNames = colnames(haplo), addEff = c(1, 2, 0.5, -1),
                   domEff = c(0.2, 0.1, 0, 0.3), epiPairs = pairs,
                   epiEff = c(1, -2))
    pop = newPop(founderPop, simParam = SP)
    gp = genParam(pop, simParam = SP)
    expect_equal(unname(gp$gv[, 1]),
                 unname(gp$mu[1] + gp$bv[, 1] + gp$dd[, 1] + gp$aa[, 1]),
                 tolerance = 1e-10)
    expect_true(all(is.finite(gp$alpha[[1]])))
    gp
  }
  forward = fit(rbind(c("m1", "m2"), c("m3", "m4")))
  reverse = fit(rbind(c("m2", "m1"), c("m4", "m3")))
  expect_equal(forward$alpha, reverse$alpha, tolerance = 1e-10)
  expect_equal(forward$bv, reverse$bv, tolerance = 1e-10)
  expect_equal(forward$dd, reverse$dd, tolerance = 1e-10)
})

test_that("THREADS genParam for epistatic traits does not depend on the thread count", {
  skip_if(getNumThreads() < 2)
  values = lapply(c(1L, 2L), function(nThreads) {
    set.seed(313)
    founderPop = runMacs(nInd = 200, nChr = 2, segSites = 100)
    SP = SimParam$new(founderPop)
    SP$nThreads = nThreads
    SP$addTraitAE(50, relAA = 0.5)
    SP$addTraitADE(50, meanDD = 0.3, relAA = 0.5)
    pop = newPop(founderPop, simParam = SP)
    genParam(pop, simParam = SP)[c("bv", "dd", "aa", "alpha", "varA", "varD", "varAA")]
  })
  expect_identical(values[[1]], values[[2]])
})

test_that("GROUP traits sharing loci keep the statistics they have alone", {
  # genParam works out correlated traits together. Their statistics for
  # each locus come from the same genotype counts either way, so they must
  # be identical to working out each trait on its own.
  set.seed(304)
  founderPop = quickHaplo(nInd = 20, nChr = 2, segSites = 200)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitAD(50, mean = rep(0, 4), var = rep(1, 4), meanDD = rep(0.2, 4),
                corA = diag(4) * 0.5 + 0.5)
  pop = randCross(newPop(founderPop, simParam = SP), 100, simParam = SP)
  grouped = AlphaSimR:::calcGenParamMulti(SP$traits, pop, 1L)
  for (i in seq_len(SP$nTraits)) {
    alone = AlphaSimR:::calcGenParam(SP$traits[[i]], pop, 1L)
    for (k in c("genicVarA", "genicVarD", "genicVarA2", "genicVarD2",
                "mu", "mu_HWE", "alpha", "alpha_HW")) {
      expect_identical(grouped[[i]][[k]], alone[[k]])
    }
    for (k in c("gv", "bv", "dd", "gv_a", "gv_d")) {
      expect_equal(grouped[[i]][[k]], alone[[k]])
    }
  }
})

test_that("THREADS genetic parameters do not depend on the thread count", {
  skip_if(getNumThreads() < 2)
  set.seed(305)
  founderPop = quickHaplo(nInd = 20, nChr = 2, segSites = 200)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitAD(200, mean = c(0, 0), var = c(1, 1), meanDD = c(0.2, 0.1),
                corA = matrix(c(1, 0.5, 0.5, 1), 2))
  SP$addTraitA(20)
  SP$addTraitADE(20, meanDD = 0.2, relAA = 0.5)
  pop = randCross(newPop(founderPop, simParam = SP), 300, simParam = SP)
  expect_identical(genParam(pop, simParam = SP, nThreads = 1L),
                   genParam(pop, simParam = SP, nThreads = 2L))
})
