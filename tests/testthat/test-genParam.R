context("genParam")

# See also test-addTrait.R & test-importData.R for some other quan. gen. tests

# Tests below are against:
# * Falconer (1961) Introduction to quantitative genetics, chapter 7, example
#   with two different allele frequencies, but assuming random mating, while
#   we here test with and without inbreeding
#   https://archive.org/details/introductiontoqu0000falc/page/112
# * TODO

test_that("genParam", {
  # rm(list=ls())
  # Haplotypes (a bit more than usual so get desired allele and genotype freqs)
  haplo = matrix(data = 0, nrow = 120, ncol = 5)

  # Locus 1
  haplo[1:8, 1] = 0 # 4 individuals with 0/0
  haplo[9:16, 1] = rep(c(0, 1), times = 4) # 4 individuals with 0/1
  haplo[17:120, 1] = 1 # 52 individuals with 1/1
  # Locus 2
  haplo[1:48, 2] = 0 # 24 individuals with 0/0
  haplo[49:120, 2] = 1 # 96 individuals with 1/1
  # Locus 3
  haplo[1:56, 3] = 0 # 28 individuals with 0/0
  haplo[57:64, 3] = rep(c(0, 1), times = 4) # 4 individuals with 0/1
  haplo[65:120, 3] = 1 # 88 individuals with 1/1
  # Locus 4
  haplo[1:30, 4] = 0 # 15 individuals with 0/0
  haplo[31:90, 4] = rep(c(0, 1), times = 30) # 30 individuals with 0/1
  haplo[91:120, 4] = 1 # 15 individuals with 1/1

  nInd = nrow(haplo) / 2
  nLoc = ncol(haplo)
  colnames(haplo) = letters[1:nLoc]

  p = colMeans(haplo)
  #   a   b   c   d   e
  # 0.9 0.6 0.5 0.5 TODO
  q = 1 - p
  #   a   b   c   d   e
  # 0.1 0.4 0.5 0.5 TODO
  P_HW = p^2
  #    a    b    c    d    e
  # 0.81 0.36 0.25 0.25 TODO
  H_HW = 2 * p * q
  #    a    b    c    d    e
  # 0.18 0.48 0.50 0.50 TODO
  Q_HW = q^2
  #    a    b    c    d    e
  # 0.01 0.16 0.25 0.25 TODO

  geno <- do.call(
    rbind,
    lapply(seq(1, nrow(haplo), by = 2), function(i) haplo[i, ] + haplo[i + 1, ])
  )
  P = apply(geno, 2, function(x) sum(x == 2) / length(x))
  #    a    b    c    d    e
  # 0.87 0.60 0.47 0.25 TODO
  H = apply(geno, 2, function(x) sum(x == 1) / length(x))
  #    a    b    c    d    e
  # 0.07 0.00 0.07 0.50 TODO
  Q = apply(geno, 2, function(x) sum(x == 0) / length(x))
  #    a    b    c    d    e
  # 0.07 0.40 0.47 0.25 TODO

  # F = (ExpHet - ObsHet) / ExpHet
  # Falconer (1961) Introduction to quantitative genetics
  # https://archive.org/details/introductiontoqu0000falc/page/66
  # page 62, formula 3.15
  F = (H_HW - H) / H_HW
  #    a    b    c    d    e
  # 0.63 1.00 0.87 0.00 TODO

  genMap = data.frame(
    markerName = letters[1:nLoc],
    chromosome = c(1, 1, 1, 1, 1),
    position = c(0, 0.2, 0.4, 0.6, 0.8)
  )

  ped = data.frame(
    id = as.character(1:nInd),
    mother = rep(0, nInd),
    father = rep(0, nInd)
  )

  founderPop = importHaplo(
    haplo = haplo,
    genMap = genMap,
    ploidy = 2L,
    ped = ped
  )

  SP = SimParam$new(founderPop = founderPop)
  SP$nThreads = 1L

  a = c(4, 4, 0, 0, 0)
  d = c(2, 2, 0, 0, 0)
  SP$importTrait(
    markerNames = "a",
    addEff = a[1],
    domEff = d[1],
    intercept = 10,
    name = "Falconer7.1_q=0.1"
  )
  SP$importTrait(
    markerNames = "b",
    addEff = a[2],
    domEff = d[2],
    intercept = 10,
    name = "Falconer7.1_q=0.4"
  )
  # SP$traits[[1]]
  # SP$traits[[2]]
  pop = newPop(founderPop, simParam = SP)
  gp = genParam(pop, simParam = SP)

  # Genetic value
  # Falconer (1961) https://archive.org/details/introductiontoqu0000falc/page/114
  expect_equal(gp$gv[c(1, 5, 9), 1], c(6, 12, 14), tolerance = 1e-6)
  expect_equal(gp$gv[c(1, 25), 2], c(6, 14), tolerance = 1e-6)

  # Origin
  # Falconer (1961) https://archive.org/details/introductiontoqu0000falc/page/114
  expect_equal(unname(gp$gv_mu[1]), 10, tolerance = 1e-6)
  expect_equal(unname(gp$gv_mu[2]), 10, tolerance = 1e-6)

  # Genetic value - part from additive effect (a) only
  # Falconer (1961) https://archive.org/details/introductiontoqu0000falc/page/114
  expect_equal(gp$gv_a[c(1, 5, 9), 1], c(-4, 0, 4), tolerance = 1e-6)
  expect_equal(gp$gv_a[c(1, 25), 2], c(-4, 4), tolerance = 1e-6)

  # Genetic value - part from dominance effect (d) only
  # Falconer (1961) https://archive.org/details/introductiontoqu0000falc/page/114
  expect_equal(gp$gv_d[c(1, 5, 9), 1], c(0, 2, 0), tolerance = 1e-6)
  expect_equal(gp$gv_d[c(1, 25), 2], c(0, 0), tolerance = 1e-6)

  # Genetic value - part from epistatic effect (aa) only
  expect_equal(gp$gv_aa[c(1, 5, 9), 1], c(0, 0, 0), tolerance = 1e-6)
  expect_equal(gp$gv_aa[c(1, 25), 2], c(0, 0), tolerance = 1e-6)

  # Genetic value - part from imprinting effect (i) only
  # TODO

  # Genetic value - part from genotype-by-environment effect (g) only
  # TODO

  # Trait mean (under random mating)
  # Origin + M in Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/115
  expect_equal(unname(gp$mu_HW[1]), 13.56, tolerance = 1e-6)
  expect_equal(unname(gp$mu_HW[2]), 11.76, tolerance = 1e-6)

  # Trait mean (actual)
  expect_equal(unname(gp$mu[1]), mean(gp$gv[, 1]), tolerance = 1e-6)
  expect_equal(unname(gp$mu[2]), mean(gp$gv[, 2]), tolerance = 1e-6)

  # Allele substitution effect (under random mating)
  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/120
  expect_equal(gp$alpha_HW[[1]][1, 1], 2.4, tolerance = 1e-6)
  expect_equal(gp$alpha_HW[[2]][1, 1], 3.6, tolerance = 1e-6)

  # Allele substitution effect (actual)
  # ... via regression
  # Fisher (1919) The Correlation between Relatives on the Supposition of Mendelian Inheritance
  # https://doi.org/10.1017/S0080456800012163
  expect_equal(
    gp$alpha[[1]][1, 1],
    cov(gp$gv[, 1], geno[, 1]) / var(geno[, 1]),
    tolerance = 1e-6
  )
  expect_equal(
    gp$alpha[[2]][1, 1],
    cov(gp$gv[, 2], geno[, 2]) / var(geno[, 2]),
    tolerance = 1e-6
  )
  # ... via formula
  # Falconer (1985) A note on Fisher’s average effect and average excess
  # https://doi.org/10.1017/S0016672300022825, page 341, formula 11
  alpha = a + d * (q - p) * ((1 - F) / (1 + F))
  expect_equal(gp$alpha[[1]][1, 1], unname(alpha[1]), tolerance = 1e-6)
  expect_equal(gp$alpha[[2]][1, 1], unname(alpha[2]), tolerance = 1e-6)

  # Allele substitution effect (actual - under random mating)
  # Antonios et al. (2025): Genetic inbreeding load and its individual prediction for milk yield in French dairy sheep
  # https://doi.org/10.1186/s12711-024-00945-z, page 2
  # TODO: add gp$alpha_F or gp$alpha_diff to genParam() or just skip it?
  alpha_F = -2 * (F / (1 + F)) * d * (q - p)
  expect_equal(
    gp$alpha[[1]][1, 1] - gp$alpha_HW[[1]][1, 1],
    unname(alpha_F[1]),
    tolerance = 1e-6
  )
  expect_equal(
    gp$alpha[[2]][1, 1] - gp$alpha_HW[[2]][1, 1],
    unname(alpha_F[2]),
    tolerance = 1e-6
  )

  # Breeding value (under random mating)
  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/121
  # TODO: add gp$bv_HW to genParam() or just skip it?
  gp$bv_HW = gp$bv
  gp$bv_HW[] = 0
  gp$bv_HW[, 1] = (geno[, 1] - 2 * p[1]) * gp$alpha_HW[[1]][1, 1]
  myBv_HW1 = unname(
    c(-2 * p[1], q[1] - p[1], 2 * q[1]) * gp$alpha_HW[[1]][1, 1]
  )
  expect_equal(myBv_HW1, c(-4.32, -1.92, 0.48), tolerance = 1e-6)
  expect_equal(gp$bv_HW[c(1, 5, 9), 1], c(-4.32, -1.92, 0.48), tolerance = 1e-6)
  meanBv = sum(myBv_HW1 * c(Q_HW[1], H_HW[1], P_HW[1]))
  expect_equal(meanBv, 0, tolerance = 1e-6)
  expect_equal(mean(gp$bv_HW[, 1]), 0, tolerance = 1e-6)

  gp$bv_HW[, 2] = (geno[, 2] - 2 * p[2]) * gp$alpha_HW[[2]][1, 1]
  myBv_HW2 = unname(
    c(-2 * p[2], q[2] - p[2], 2 * q[2]) * gp$alpha_HW[[2]][1, 1]
  )
  expect_equal(myBv_HW2, c(-4.32, -0.72, 2.88), tolerance = 1e-6)
  expect_equal(gp$bv_HW[c(1, 25), 2], c(-4.32, 2.88), tolerance = 1e-6)
  meanBv = sum(myBv_HW2 * c(Q_HW[2], H_HW[2], P_HW[2]))
  expect_equal(meanBv, 0, tolerance = 1e-6)
  expect_equal(mean(gp$bv_HW[, 2]), 0, tolerance = 1e-6)

  # Breeding value (actual)
  myBv = unname(c(-2 * p[1], q[1] - p[1], 2 * q[1]) * gp$alpha[[1]][1, 1])
  expect_equal(gp$bv[c(1, 5, 9), 1], myBv, tolerance = 1e-6)
  meanBv = sum(myBv * c(Q[1], H[1], P[1]))
  expect_equal(meanBv, 0, tolerance = 1e-6)
  expect_equal(mean(gp$bv[, 1]), 0, tolerance = 1e-6)

  myBv = unname(c(-2 * p[2], q[2] - p[2], 2 * q[2]) * gp$alpha[[2]][1, 1])
  expect_equal(gp$bv[c(1, 25), 2], myBv[c(1, 3)], tolerance = 1e-6)
  meanBv = sum(myBv * c(Q[2], H[2], P[2]))
  expect_equal(meanBv, 0, tolerance = 1e-6)
  expect_equal(mean(gp$bv[, 2]), 0, tolerance = 1e-6)

  # Inbreeding depression load
  # Antonios et al. (2025): Genetic inbreeding load and its individual prediction for milk yield in French dairy sheep
  # https://doi.org/10.1186/s12711-024-00945-z, page 2 and 3
  # TODO: add gp$idl to genParam() or just skip it?
  gp$idl = gp$bv - gp$bv_HW
  myIdl = unname(c(-2 * p[1], q[1] - p[1], 2 * q[1]) * alpha_F[1])
  expect_equal(gp$idl[c(1, 5, 9), 1], myIdl, tolerance = 1e-6)
  myIdl = unname(c(-2 * p[2], q[2] - p[2], 2 * q[2]) * alpha_F[2])
  expect_equal(gp$idl[c(1, 25), 2], myIdl[c(1, 3)], tolerance = 1e-6)

  # Dominance deviation (random mating)
  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/123
  gp$dd_HW = gp$dd
  gp$dd_HW[] = 0
  gp$dd_HW[, 1] = gp$gv[, 1] - gp$mu_HW[1] - gp$bv_HW[, 1]
  # gp$gv[c(1, 5, 9), 1]
  # 6 12 14
  # gp$gv[c(1, 5, 9), 1] - gp$mu_HW[1]
  # -7.56 -1.56  0.44
  # gp$bv_HW[c(1, 5, 9), 1]
  # -4.32 -1.92  0.48
  # gp$dd_HW[c(1, 5, 9), 1]
  # -3.24  0.36 -0.04
  # This is c(-2*p[1]^2*d[1], 2*p[1]*q[1]*d[1], -2*q[1]^2*d[1])
  expect_equal(gp$dd_HW[c(1, 5, 9), 1], c(-3.24, 0.36, -0.04), tolerance = 1e-6)
  gp$dd_HW[, 2] = gp$gv[, 2] - gp$mu_HW[2] - gp$bv_HW[, 2]
  # gp$gv[c(1, 25), 2]
  # 6 14
  # gp$gv[c(1, 25), 2] - gp$mu_HW[2]
  # -5.76  2.24
  # gp$bv_HW[c(1, 25), 2]
  # -4.32 -1.92  0.48
  # gp$dd_HW[c(1, 25), 2]
  # -1.44 -0.64
  # This is c(-2*p[1]^2*d[1], -2*q[1]^2*d[1])
  expect_equal(gp$dd_HW[c(1, 25), 2], c(-1.44, -0.64), tolerance = 1e-6)

  # Dominance deviation (actual)
  # gp$gv[c(1, 5, 9), 1]
  # 6 12 14
  # gp$gv[c(1, 5, 9), 1] - gp$mu[1]
  # -7.3333333 -1.3333333  0.6666667
  # gp$bv[c(1, 5, 9), 1]
  # -6.5454545 -2.9090909  0.7272727
  # gp$dd[c(1, 5, 9), 1]
  # -0.78787879  1.57575758 -0.06060606
  myDd = gp$gv[c(1, 5, 9), 1] - gp$mu[1] - gp$bv[c(1, 5, 9), 1]
  # This is not the same as c(-2*p[1]^2*d[1], 2*p[1]*q[1]*d[1], -2*q[1]^2*d[1])
  # due to inbreeding!
  expect_equal(gp$dd[c(1, 5, 9), 1], myDd, tolerance = 1e-6)
  # gp$gv[c(1, 25), 2]
  # 6 14
  # gp$gv[c(1, 25), 2] - gp$mu[2]
  # -4.8  3.2
  # gp$bv[c(1, 25), 2]
  # -4.8  3.2
  # gp$dd[c(1, 25), 2]
  # 0 0
  myDd = gp$gv[c(1, 25), 2] - gp$mu[2] - gp$bv[c(1, 25), 2]
  # This is not the same as c(-2*p[1]^2*d[1], -2*q[1]^2*d[1])
  # due to inbreeding!
  expect_equal(gp$dd[c(1, 25), 2], myDd, tolerance = 1e-6)

  # TODO: aa
  # TODO: VarA
  # TODO: VarD
  # TODO: varAA
  # TODO: varG
  # TODO: genicVarA
  # TODO: genicVarD
  # TODO: genicVarAA
  # TODO: genicVarG
  # TODO: covA_HW
  # TODO: covD_HW
  # TODO: covAA_HW
  # TODO: covA_L
  # TODO: covD_L
  # TODO: covAA_L
  # TODO: covAA_L
  # TODO: covAD_L
  # TODO: covAAA_L
  # TODO: covDAA_L
  # TODO: covG_L
})
