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

expectGenParamMatchesRef = function(pop, simParam) {
  gp = genParam(pop, simParam = simParam)
  for (i in seq_len(simParam$nTraits)) {
    ref = refBvDd(pop, simParam, i)
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
