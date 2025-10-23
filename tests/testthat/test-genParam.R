context("genParam")

# See also test-addTrait.R & test-importData.R for some other quan. gen. tests

# Tests below are against:
# * Falconer (1961) Introduction to quantitative genetics, chapter 7, example
#   with two different allele frequencies, but assuming Hardy-Weinberg equilibrium
#   (=random mating, no inbreeding), while we here test with and without this
#   assumption inbreeding
#   https://archive.org/details/introductiontoqu0000falc/page/112
# * TODO

test_that("genParam", {
  # rm(list=ls())
  # Haplotypes (a bit more than usual so get desired allele and genotype freqs)
  haplo = matrix(data = 0, nrow = 200, ncol = 3)

  # Locus 1 - in Hardy-Weinberg equilibrium (=random mating, no inbreeding)
  # Aiming for q = 0.1 as in Falconer (1996) example 7.1
  haplo[1:2, 1] = 0 # 1 individual with 0/0
  haplo[3:38, 1] = rep(c(0, 1), times = 18) # 18 individuals with 0/1
  haplo[39:200, 1] = 1 # 81 individuals with 1/1

  # Locus 2 - in Hardy-Weinberg equilibrium (=random mating, no inbreeding)
  # Aiming for q = 0.4 as in Falconer (1996) example 7.1
  haplo[1:32, 2] = 0 # 16 individuals with 0/0
  haplo[33:128, 2] = rep(c(0, 1), times = 48) # 48 individuals with 0/1
  haplo[129:200, 2] = 1 # 36 individuals with 1/1

  # Locus 3 - not in Hardy-Weinberg equilibrium (=non-random mating, inbreeding)
  # Aiming for q = 0.4 as in Falconer (1996) example 7.1
  haplo[1:68, 3] = 0 # 34 individuals with 0/0
  haplo[69:92, 3] = rep(c(0, 1), times = 12) # 12 individuals with 0/1
  haplo[93:200, 3] = 1 # 54 individuals with 1/1

  # TODO Locus 3 - TODO
  # Aiming for q = 0.4 as in Falconer (1996) example 7.1
  # haplo[1:48, 2] = 0 # 24 individuals with 0/0
  # haplo[49:120, 2] = 1 # 96 individuals with 1/1

  # TODO haplo, locus 3 - TODO
  # haplo[1:30, 3] = 0 # 15 individuals with 0/0
  # haplo[31:90, 3] = rep(c(0, 1), times = 30) # 30 individuals with 0/1
  # haplo[91:120, 3] = 1 # 15 individuals with 1/1

  # TODO haplo, locus TODO
  # haplo[1:56, 3] = 0 # 28 individuals with 0/0
  # haplo[57:64, 3] = rep(c(0, 1), times = 4) # 4 individuals with 0/1
  # haplo[65:120, 3] = 1 # 88 individuals with 1/1

  nInd = nrow(haplo) / 2
  nLoc = ncol(haplo)
  colnames(haplo) = letters[1:nLoc]

  p = colMeans(haplo)
  #   a   b   c   d   e
  # 0.9 0.6 0.6 TODO TODO
  q = 1 - p
  #   a   b   c   d   e
  # 0.1 0.4 0.4 TODO TODO
  P_HW = p^2
  #    a    b    c    d    e
  # 0.81 0.36 0.36 TODO TODO
  H_HW = 2 * p * q
  #    a    b    c    d    e
  # 0.18 0.48 0.48 TODO TODO
  Q_HW = q^2
  #    a    b    c    d    e
  # 0.01 0.16 0.16 TODO TODO

  geno <- do.call(
    rbind,
    lapply(seq(1, nrow(haplo), by = 2), function(i) haplo[i, ] + haplo[i + 1, ])
  )
  P = apply(geno, 2, function(x) sum(x == 2) / length(x))
  #    a    b    c    d    e
  # 0.81 0.36 0.54 TODO TODO
  # P * nInd
  #  a  b  c
  # 81 36 54
  H = apply(geno, 2, function(x) sum(x == 1) / length(x))
  #    a    b    c    d    e
  # 0.18 0.48 0.12 TODO TODO
  # H * nInd
  #  a  b  c
  # 18 48 12
  Q = apply(geno, 2, function(x) sum(x == 0) / length(x))
  #    a    b    c    d    e
  # 0.01 0.16 0.34 TODO TODO
  # Q * nInd
  #  a  b  c
  #  1 16 34

  # F = (ExpHet - ObsHet) / ExpHet
  # Falconer (1961) Introduction to quantitative genetics
  # https://archive.org/details/introductiontoqu0000falc/page/66
  # page 62, formula 3.15
  F = (H_HW - H) / H_HW
  # a b    c    d    e
  # 0 0 0.75 TODO TODO

  genMap = data.frame(
    markerName = letters[1:nLoc],
    chromosome = c(1, 1, 1),
    position = c(0, 0.2, 0.4)
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

  a = c(4, 4, 4)
  d = c(2, 2, 2)
  SP$importTrait(
    markerNames = "a",
    addEff = a[1],
    domEff = d[1],
    intercept = 10,
    name = "Falconer7.1_q=0.1_HWE"
  )
  SP$importTrait(
    markerNames = "b",
    addEff = a[2],
    domEff = d[2],
    intercept = 10,
    name = "Falconer7.1_q=0.4_HWE"
  )
  SP$importTrait(
    markerNames = "c",
    addEff = a[3],
    domEff = d[3],
    intercept = 10,
    name = "Falconer7.1_q=0.4_not_HWE"
  )
  # SP$traits[[1]]
  # SP$traits[[2]]
  # SP$traits[[3]]
  pop = newPop(founderPop, simParam = SP)
  gp = genParam(pop, simParam = SP)

  # ---- Genetic value ----

  genoLoc1 = c(1, 2, 20)
  genoLoc2 = c(1, 17, 65)
  genoLoc3 = c(1, 35, 47)

  # Falconer (1961) https://archive.org/details/introductiontoqu0000falc/page/114
  expect_equal(gp$gv[genoLoc1, 1], c(6, 12, 14), tolerance = 1e-6)
  expect_equal(gp$gv[genoLoc2, 2], c(6, 12, 14), tolerance = 1e-6)
  expect_equal(gp$gv[genoLoc3, 3], c(6, 12, 14), tolerance = 1e-6)

  # ---- Origin ----

  # Falconer (1961) https://archive.org/details/introductiontoqu0000falc/page/114
  expect_equal(unname(gp$gv_mu[1]), 10, tolerance = 1e-6)
  expect_equal(unname(gp$gv_mu[2]), 10, tolerance = 1e-6)
  expect_equal(unname(gp$gv_mu[3]), 10, tolerance = 1e-6)

  # ---- Genetic value - part from additive effect (a) only ----

  # Falconer (1961) https://archive.org/details/introductiontoqu0000falc/page/114
  expect_equal(gp$gv_a[genoLoc1, 1], c(-4, 0, 4), tolerance = 1e-6)
  expect_equal(gp$gv_a[genoLoc2, 2], c(-4, 0, 4), tolerance = 1e-6)
  expect_equal(gp$gv_a[genoLoc3, 3], c(-4, 0, 4), tolerance = 1e-6)

  # ---- Genetic value - part from dominance effect (d) only ----

  # Falconer (1961) https://archive.org/details/introductiontoqu0000falc/page/114
  expect_equal(gp$gv_d[genoLoc1, 1], c(0, 2, 0), tolerance = 1e-6)
  expect_equal(gp$gv_d[genoLoc2, 2], c(0, 2, 0), tolerance = 1e-6)
  expect_equal(gp$gv_d[genoLoc3, 3], c(0, 2, 0), tolerance = 1e-6)

  # ---- Genetic value - part from additive-by-additive effect (aa) only ----

  expect_equal(gp$gv_aa[genoLoc1, 1], c(0, 0, 0), tolerance = 1e-6)
  expect_equal(gp$gv_aa[genoLoc2, 2], c(0, 0, 0), tolerance = 1e-6)
  expect_equal(gp$gv_aa[genoLoc3, 2], c(0, 0, 0), tolerance = 1e-6)

  # ---- Genetic value - part from imprinting effect (i) only ----

  # TODO

  # ---- Genetic value - part from genotype-by-environment effect (g) only ----

  # TODO

  # ---- Trait mean (under Hardy-Weinberg equilibrium) ----

  # Origin + M in Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/115
  expect_equal(unname(gp$mu_HW[1]), 13.56, tolerance = 1e-6)
  expect_equal(unname(gp$mu_HW[2]), 11.76, tolerance = 1e-6)
  expect_equal(unname(gp$mu_HW[3]), 11.76, tolerance = 1e-6)

  # ---- Trait mean (actual) ----

  expect_equal(unname(gp$mu[1]), mean(gp$gv[, 1]), tolerance = 1e-6) # 13.56
  expect_equal(gp$mu[1], gp$mu_HW[1], tolerance = 1e-6) # 13.56
  expect_equal(unname(gp$mu[2]), mean(gp$gv[, 2]), tolerance = 1e-6) # 11.76
  expect_equal(gp$mu[2], gp$mu_HW[2], tolerance = 1e-6) # 11.76
  expect_equal(unname(gp$mu[3]), mean(gp$gv[, 3]), tolerance = 1e-6) # 11.04

  # ---- Allele substitution effect (under Hardy-Weinberg equilibrium) ----

  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/120
  expect_equal(gp$alpha_HW[[1]][1, 1], 2.4, tolerance = 1e-6)
  expect_equal(gp$alpha_HW[[2]][1, 1], 3.6, tolerance = 1e-6)
  expect_equal(gp$alpha_HW[[3]][1, 1], 3.6, tolerance = 1e-6)

  # ---- Allele substitution effect (actual) ----

  # ... via regression
  # Fisher (1919) The Correlation between Relatives on the Supposition of Mendelian Inheritance
  # https://doi.org/10.1017/S0080456800012163
  expect_equal(
    gp$alpha[[1]][1, 1],
    cov(gp$gv[, 1], geno[, 1]) / var(geno[, 1]),
    tolerance = 1e-6
  ) # 2.4
  expect_equal(
    gp$alpha[[2]][1, 1],
    cov(gp$gv[, 2], geno[, 2]) / var(geno[, 2]),
    tolerance = 1e-6
  ) # 3.6
  expect_equal(
    gp$alpha[[3]][1, 1],
    cov(gp$gv[, 3], geno[, 3]) / var(geno[, 3]),
    tolerance = 1e-6
  ) # 3.942857
  # ... via formula
  # Falconer (1985) A note on Fisher’s average effect and average excess
  # https://doi.org/10.1017/S0016672300022825, page 341, formula 11
  alpha = a + d * (q - p) * ((1 - F) / (1 + F))
  expect_equal(gp$alpha[[1]][1, 1], unname(alpha[1]), tolerance = 1e-6) # 2.4
  expect_equal(gp$alpha[[2]][1, 1], unname(alpha[2]), tolerance = 1e-6) # 3.6
  expect_equal(gp$alpha[[3]][1, 1], unname(alpha[3]), tolerance = 1e-6) # 3.942857

  # ---- Allele substitution effect (actual - under Hardy-Weinberg equilibrium) ----

  # TODO: add gp$alpha_F or gp$alpha_diff to genParam() or just skip it?
  # Antonios et al. (2025): Genetic inbreeding load and its individual prediction for milk yield in French dairy sheep
  # https://doi.org/10.1186/s12711-024-00945-z, page 2

  alpha_F = -2 * (F / (1 + F)) * d * (q - p)
  expect_equal(
    gp$alpha[[1]][1, 1] - gp$alpha_HW[[1]][1, 1],
    unname(alpha_F[1]),
    tolerance = 1e-6
  ) # 0
  expect_equal(
    gp$alpha[[2]][1, 1] - gp$alpha_HW[[2]][1, 1],
    unname(alpha_F[2]),
    tolerance = 1e-6
  ) # 0
  expect_equal(
    gp$alpha[[3]][1, 1] - gp$alpha_HW[[3]][1, 1],
    unname(alpha_F[3]),
    tolerance = 1e-6
  ) # 0.342857 = 3.942857 - 3.6

  # ---- Breeding value (under Hardy-Weinberg equilibrium) ----

  # TODO: Name this as Additive genetic value?
  # TODO: add gp$bv_HW to genParam() or just skip it?

  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/121
  gp$bv_HW = gp$bv
  gp$bv_HW[] = 0
  gp$bv_HW[, 1] = (geno[, 1] - 2 * p[1]) * gp$alpha_HW[[1]][1, 1]
  myBv_HW1 = unname(
    c(-2 * p[1], q[1] - p[1], 2 * q[1]) * gp$alpha_HW[[1]][1, 1]
  )
  expect_equal(myBv_HW1, c(-4.32, -1.92, 0.48), tolerance = 1e-6)
  expect_equal(
    gp$bv_HW[genoLoc1, 1],
    c(-4.32, -1.92, 0.48),
    tolerance = 1e-6
  )
  meanBv = sum(myBv_HW1 * c(Q_HW[1], H_HW[1], P_HW[1]))
  expect_equal(meanBv, 0, tolerance = 1e-6)
  expect_equal(mean(gp$bv_HW[, 1]), 0, tolerance = 1e-6)

  gp$bv_HW[, 2] = (geno[, 2] - 2 * p[2]) * gp$alpha_HW[[2]][1, 1]
  myBv_HW2 = unname(
    c(-2 * p[2], q[2] - p[2], 2 * q[2]) * gp$alpha_HW[[2]][1, 1]
  )
  expect_equal(myBv_HW2, c(-4.32, -0.72, 2.88), tolerance = 1e-6)
  expect_equal(
    gp$bv_HW[genoLoc2, 2],
    c(-4.32, -0.72, 2.88),
    tolerance = 1e-6
  )
  meanBv = sum(myBv_HW2 * c(Q_HW[2], H_HW[2], P_HW[2]))
  expect_equal(meanBv, 0, tolerance = 1e-6)
  expect_equal(mean(gp$bv_HW[, 2]), 0, tolerance = 1e-6)

  gp$bv_HW[, 3] = (geno[, 3] - 2 * p[3]) * gp$alpha_HW[[3]][1, 1]
  myBv_HW3 = unname(
    c(-2 * p[3], q[3] - p[3], 2 * q[3]) * gp$alpha_HW[[3]][1, 1]
  )
  expect_equal(myBv_HW3, c(-4.32, -0.72, 2.88), tolerance = 1e-6)
  expect_equal(
    gp$bv_HW[genoLoc3, 3],
    c(-4.32, -0.72, 2.88),
    tolerance = 1e-6
  )
  meanBv = sum(myBv_HW3 * c(Q_HW[3], H_HW[3], P_HW[3]))
  expect_equal(meanBv, 0, tolerance = 1e-6)
  expect_equal(mean(gp$bv_HW[, 3]), 0, tolerance = 1e-6)

  # ---- Breeding value (actual) ----

  # TODO: Name this as Additive genetic value?
  myBv = unname(c(-2 * p[1], q[1] - p[1], 2 * q[1]) * gp$alpha[[1]][1, 1])
  expect_equal(gp$bv[genoLoc1, 1], myBv, tolerance = 1e-6)
  meanBv = sum(myBv * c(Q[1], H[1], P[1]))
  expect_equal(meanBv, 0, tolerance = 1e-6)
  expect_equal(mean(gp$bv[, 1]), 0, tolerance = 1e-6)

  myBv = unname(c(-2 * p[2], q[2] - p[2], 2 * q[2]) * gp$alpha[[2]][1, 1])
  expect_equal(gp$bv[genoLoc2, 2], myBv[c(1, 2, 3)], tolerance = 1e-6)
  meanBv = sum(myBv * c(Q[2], H[2], P[2]))
  expect_equal(meanBv, 0, tolerance = 1e-6)
  expect_equal(mean(gp$bv[, 2]), 0, tolerance = 1e-6)

  myBv = unname(c(-2 * p[3], q[3] - p[3], 2 * q[3]) * gp$alpha[[3]][1, 1])
  expect_equal(gp$bv[genoLoc3, 3], myBv[c(1, 2, 3)], tolerance = 1e-6)
  meanBv = sum(myBv * c(Q[3], H[3], P[3]))
  expect_equal(meanBv, 0, tolerance = 1e-6)
  expect_equal(mean(gp$bv[, 3]), 0, tolerance = 1e-6)

  # ---- Inbreeding depression load ----

  # TODO: add gp$idl to genParam() or just skip it?

  # Antonios et al. (2025): Genetic inbreeding load and its individual prediction for milk yield in French dairy sheep
  # https://doi.org/10.1186/s12711-024-00945-z, page 2 and 3

  gp$idl = gp$bv - gp$bv_HW
  myIdl = unname(c(-2 * p[1], q[1] - p[1], 2 * q[1]) * alpha_F[1])
  expect_equal(gp$idl[genoLoc1, 1], myIdl, tolerance = 1e-6)
  myIdl = unname(c(-2 * p[2], q[2] - p[2], 2 * q[2]) * alpha_F[2])
  expect_equal(gp$idl[genoLoc2, 2], myIdl, tolerance = 1e-6)
  myIdl = unname(c(-2 * p[3], q[3] - p[3], 2 * q[3]) * alpha_F[3])
  expect_equal(gp$idl[genoLoc3, 3], myIdl, tolerance = 1e-6)

  # ---- Dominance deviation (under Hardy-Weinberg equilibrium) ----

  # TODO: add gp$bv_HW to genParam() or just skip it?

  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/123
  gp$dd_HW = gp$dd
  gp$dd_HW[] = 0

  gp$dd_HW[, 1] = gp$gv[, 1] - gp$mu_HW[1] - gp$bv_HW[, 1]
  # gp$gv[genoLoc1, 1]
  # 6 12 14
  # gp$gv[genoLoc1, 1] - gp$mu_HW[1]
  # -7.56 -1.56  0.44
  # gp$bv_HW[genoLoc1, 1]
  # -4.32 -1.92  0.48
  # gp$dd_HW[genoLoc1, 1]
  # -3.24  0.36 -0.04
  # This is c(-2*p[1]^2*d[1], 2*p[1]*q[1]*d[1], -2*q[1]^2*d[1])
  expect_equal(gp$dd_HW[genoLoc1, 1], c(-3.24, 0.36, -0.04), tolerance = 1e-6)

  gp$dd_HW[, 2] = gp$gv[, 2] - gp$mu_HW[2] - gp$bv_HW[, 2]
  # gp$gv[genoLoc2, 2]
  # 6 12 14
  # gp$gv[genoLoc2, 2] - gp$mu_HW[2]
  # -5.76 0.24 2.24
  # gp$bv_HW[genoLoc2, 2]
  # -4.32 -0.72  2.88
  # gp$dd_HW[genoLoc2, 2]
  # -1.44 0.96 -0.64
  # This is c(-2*p[2]^2*d[2], 2*p[2]*q[2]*d[2], -2*q[2]^2*d[2])
  expect_equal(gp$dd_HW[genoLoc2, 2], c(-1.44, 0.96, -0.64), tolerance = 1e-6)

  gp$dd_HW[, 3] = gp$gv[, 3] - gp$mu_HW[3] - gp$bv_HW[, 3]
  # gp$gv[genoLoc3, 3]
  # 6 12 14
  # gp$gv[genoLoc3, 3] - gp$mu_HW[3]
  # -5.76 0.24 2.24
  # gp$bv_HW[genoLoc3, 3]
  # -4.32 -0.72  2.88
  # gp$dd_HW[genoLoc3, 3]
  # -1.44 0.96 -0.64
  # This is c(-2*p[3]^2*d[3], 2*p[3]*q[3]*d[3], -2*q[3]^2*d[3])
  expect_equal(gp$dd_HW[genoLoc3, 3], c(-1.44, 0.96, -0.64), tolerance = 1e-6)

  # ---- Dominance deviation (actual) ----

  # gp$gv[genoLoc1, 1]
  # 6 12 14
  # gp$gv[genoLoc1, 1] - gp$mu[1]
  # -7.56 -1.56  0.44
  # gp$bv[genoLoc1, 1]
  # -4.32 -1.92  0.48
  # gp$dd[genoLoc1, 1]
  # -3.24 0.36 -0.04
  myDd = gp$gv[genoLoc1, 1] - gp$mu[1] - gp$bv[genoLoc1, 1]
  # This is the same as c(-2*p[1]^2*d[1], 2*p[1]*q[1]*d[1], -2*q[1]^2*d[1]) due to HWE
  expect_equal(gp$dd[genoLoc1, 1], myDd, tolerance = 1e-6)
  # gp$gv[genoLoc2, 2]
  # 6 12 14
  # gp$gv[genoLoc2, 2] - gp$mu[2]
  # -5.76  0.24  2.24
  # gp$bv[genoLoc2, 2]
  # -4.32 -0.72  2.88
  # gp$dd[genoLoc2, 2]
  # -1.44  0.96 -0.64
  myDd = gp$gv[genoLoc2, 2] - gp$mu[2] - gp$bv[genoLoc2, 2]
  # This is the same as c(-2*p[2]^2*d[2], 2*p[2]*q[2]*d[2], -2*q[2]^2*d[2]) due to HWE
  expect_equal(gp$dd[genoLoc2, 2], myDd, tolerance = 1e-6)

  # gp$gv[genoLoc1, 1]
  # 6 12 14
  # gp$gv[genoLoc1, 1] - gp$mu[1]
  # -7.56 -1.56  0.44
  # gp$bv[genoLoc1, 1]
  # -4.32 -1.92  0.48
  # gp$dd[genoLoc1, 1]
  # -3.24 0.36 -0.04
  myDd = gp$gv[genoLoc1, 1] - gp$mu[1] - gp$bv[genoLoc1, 1]
  # This is the same as c(-2*p[1]^2*d[1], 2*p[1]*q[1]*d[1], -2*q[1]^2*d[1]) due to HWE
  expect_equal(gp$dd[genoLoc1, 1], myDd, tolerance = 1e-6)

  # gp$gv[genoLoc3, 3]
  # 6 12 14
  # gp$gv[genoLoc3, 3] - gp$mu[3]
  # -5.04  0.96  2.96
  # gp$bv[genoLoc3, 3]
  # -4.7314286 -0.7885714  3.1542857
  # gp$dd[genoLoc3, 3]
  # -0.3085714  1.7485714 -0.1942857
  myDd = gp$gv[genoLoc3, 3] - gp$mu[3] - gp$bv[genoLoc3, 3]
  # This is NOT the same as c(-2*p[3]^2*d[3], 2*p[3]*q[3]*d[3], -2*q[3]^2*d[3]) due to inbreeding
  expect_equal(gp$dd[genoLoc3, 3], myDd, tolerance = 1e-6)

  # ---- Imprinting deviation (under Hardy-Weinberg equilibrium) ----

  # TODO

  # ---- Imprinting deviation (actual) ----

  # TODO

  # ---- Additive-by-additive deviation (under Hardy-Weinberg equilibrium) ----

  # No AxA epistasis in loci 1-3
  expect_equal(gp$aa[genoLoc1, 1], c(0, 0, 0), tolerance = 1e-6)
  expect_equal(gp$aa[genoLoc2, 2], c(0, 0, 0), tolerance = 1e-6)
  expect_equal(gp$aa[genoLoc3, 2], c(0, 0, 0), tolerance = 1e-6)

  # ---- Additive-by-additive deviation (actual) ----

  # No AxA epistasis in loci 1-3
  expect_equal(gp$aa[genoLoc1, 1], c(0, 0, 0), tolerance = 1e-6)
  expect_equal(gp$aa[genoLoc2, 2], c(0, 0, 0), tolerance = 1e-6)
  expect_equal(gp$aa[genoLoc3, 2], c(0, 0, 0), tolerance = 1e-6)

  # ---- (Total) Genetic variance (under Hardy-Weinberg equilibrium) ----

  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/136
  myVarG = sum(
    (gp$gv[genoLoc1, 1] - gp$mu_HW[1])^2 * c(Q_HW[1], H_HW[1], P_HW[1])
  )
  expect_equal(myVarG, 1.1664, tolerance = 1e-6)
  myVarG = 2 *
    p[1] *
    q[1] *
    gp$alpha_HW[[1]][1, 1]^2 +
    (2 * p[1] * q[1] * d[1])^2
  expect_equal(unname(myVarG), 1.1664, tolerance = 1e-6)

  myVarG = sum(
    (gp$gv[genoLoc2, 2] - gp$mu_HW[2])^2 * c(Q_HW[2], H_HW[2], P_HW[2])
  )
  expect_equal(myVarG, 7.1424, tolerance = 1e-6)
  myVarG = 2 *
    p[2] *
    q[2] *
    gp$alpha_HW[[2]][1, 1]^2 +
    (2 * p[2] * q[2] * d[2])^2
  expect_equal(unname(myVarG), 7.1424, tolerance = 1e-6)

  myVarG = sum(
    (gp$gv[genoLoc3, 3] - gp$mu_HW[3])^2 * c(Q_HW[3], H_HW[3], P_HW[3])
  )
  expect_equal(myVarG, 7.1424, tolerance = 1e-6)
  myVarG = 2 *
    p[3] *
    q[3] *
    gp$alpha_HW[[3]][1, 1]^2 +
    (2 * p[3] * q[3] * d[3])^2
  expect_equal(unname(myVarG), 7.1424, tolerance = 1e-6)

  # ---- (Total) Genetic variance (actual) ----

  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/136
  myVarG = sum((gp$gv[genoLoc1, 1] - gp$mu[1])^2 * c(Q[1], H[1], P[1])) # 1.1664
  expect_equal(gp$varG[1, 1], myVarG, tolerance = 1e-6)
  expect_equal(
    gp$varG[1, 1],
    popVar(gp$gv[, 1, drop = FALSE])[1, 1],
    tolerance = 1e-6
  )

  myVarG = sum((gp$gv[genoLoc2, 2] - gp$mu[2])^2 * c(Q[2], H[2], P[2])) # 7.1424
  expect_equal(gp$varG[2, 2], myVarG, tolerance = 1e-6)
  expect_equal(
    gp$varG[2, 2],
    popVar(gp$gv[, 2, drop = FALSE])[1, 1],
    tolerance = 1e-6
  )

  myVarG = sum((gp$gv[genoLoc3, 3] - gp$mu[3])^2 * c(Q[3], H[3], P[3])) # 13.4784
  expect_equal(gp$varG[3, 3], myVarG, tolerance = 1e-6)
  expect_equal(
    gp$varG[3, 3],
    popVar(gp$gv[, 3, drop = FALSE])[1, 1],
    tolerance = 1e-6
  )

  # ---- Additive genetic variance (under Hardy-Weinberg equilibrium) ----

  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/136
  theoVarA = unname(2 * p[1] * q[1] * gp$alpha_HW[[1]][1, 1]^2)
  expect_equal(
    theoVarA,
    popVar(gp$bv_HW[, 1, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 1.0368
  expect_equal(theoVarA, gp$varA[1, 1], tolerance = 1e-6) # 1.0368
  theoVarA = unname(2 * p[2] * q[2] * gp$alpha_HW[[2]][1, 1]^2)
  expect_equal(
    theoVarA,
    popVar(gp$bv_HW[, 2, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 6.2208
  expect_equal(theoVarA, gp$varA[2, 2], tolerance = 1e-6) # 6.2208
  theoVarA = unname(2 * p[3] * q[3] * gp$alpha_HW[[3]][1, 1]^2) # 6.2208
  myVarA = sum(gp$bv_HW[genoLoc3, 3]^2 * c(Q_HW[3], H_HW[3], P_HW[3])) # 6.2208
  expect_equal(
    unname(2 * p[3] * q[3] * gp$alpha_HW[[3]][1, 1]^2), # 6.2208
    myVarA,
    # popVar(gp$bv_HW[, 3, drop = FALSE])[1, 1], # 10.8864
    # Due to deviation from HWE
    tolerance = 1e-6
  )
  # FAILS expect_equal(myVarA, # 6.2208
  #                    unname(gp$varA[3,3]), tolerance = 1e-6) # 13.05874
  # Due to deviation from HWE

  # ---- Additive genetic variance (actual) ----

  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/136
  theoVarA = unname(2 * p[1] * q[1] * gp$alpha[[1]][1, 1]^2 * (1 + F[1]))
  expect_equal(
    theoVarA,
    popVar(gp$bv[, 1, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 1.0368
  expect_equal(theoVarA, gp$varA[1, 1], tolerance = 1e-6) # 1.0368
  theoVarA = unname(2 * p[2] * q[2] * gp$alpha[[2]][1, 1]^2 * (1 + F[2]))
  expect_equal(
    theoVarA,
    popVar(gp$bv[, 2, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 6.2208
  expect_equal(theoVarA, gp$varA[2, 2], tolerance = 1e-6) # 1.0368
  myVarA = sum(gp$bv[genoLoc3, 3]^2 * c(Q[3], H[3], P[3])) # 13.05874
  expect_equal(
    # unname(2 * p[3] * q[3] * gp$alpha[[3]][1, 1]^2), # 7.462139
    # We don't have HWE, so the above does not hold
    unname(2 * p[3] * q[3] * gp$alpha[[3]][1, 1]^2 * (1 + F[3])), # 13.05874
    popVar(gp$bv[, 3, drop = FALSE])[1, 1], # 13.05874
    myVarA,
    tolerance = 1e-6
  )
  expect_equal(myVarA, gp$varA[3, 3], tolerance = 1e-6) # 13.05874

  # ---- Dominance genetic variance (under Hardy-Weinberg equilibrium) ----

  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/136
  theoVarD = unname((2 * p[1] * q[1] * d[1])^2)
  expect_equal(
    theoVarD,
    popVar(gp$dd_HW[, 1, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 0.1296
  expect_equal(theoVarD, gp$varD[1, 1], tolerance = 1e-6) # 0.1296
  theoVarD = unname((2 * p[2] * q[2] * d[2])^2)
  expect_equal(
    theoVarD,
    popVar(gp$dd_HW[, 2, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 0.9216
  expect_equal(theoVarD, gp$varD[2, 2], tolerance = 1e-6) # 0.9216
  theoVarD = unname((2 * p[3] * q[3] * d[3])^2) # 0.9216
  myVarD = sum(gp$dd_HW[genoLoc3, 3]^2 * c(Q_HW[3], H_HW[3], P_HW[3])) # 0.9216
  expect_equal(
    unname((2 * p[3] * q[3] * d[3])^2), # 0.9216
    myVarD,
    # popVar(gp$dd_HW[, 3, drop = FALSE])[1, 1], # 0.5184
    # Due to deviation from HWE
    tolerance = 1e-6
  )
  # FAILS expect_equal(myVarD, # 0.9216
  #                    unname(gp$varD[3,3]), tolerance = 1e-6) # 0.4196571
  # Due to deviation from HWE

  # ---- Dominance genetic variance (actual) ----

  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/136
  theoVarD = unname((2 * p[1] * q[1] * d[1])^2)
  expect_equal(
    theoVarD,
    popVar(gp$dd[, 1, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 0.1296
  expect_equal(theoVarD, gp$varD[1, 1], tolerance = 1e-6) # 0.1296
  theoVarD = unname((2 * p[2] * q[2] * d[2])^2)
  expect_equal(
    theoVarD,
    popVar(gp$dd[, 2, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 0.9216
  expect_equal(theoVarD, gp$varD[2, 2], tolerance = 1e-6) # 0.9216
  myVarD = sum(gp$dd[genoLoc3, 3]^2 * c(Q[3], H[3], P[3])) # 0.4196571
  expect_equal(
    # unname((2 * p[3] * q[3] * d[3])^2), # 0.9216
    # We don't have HWE, so the above does not hold
    popVar(gp$dd[, 3, drop = FALSE])[1, 1], # 0.4196571
    myVarD,
    tolerance = 1e-6
  )
  theoVarD = unname((2 * p[3] * q[3] * d[3])^2) # 0.9216
  myVarD = sum(gp$dd[genoLoc3, 3]^2 * c(Q[3], H[3], P[3])) # 0.4196571
  expect_equal(
    # unname((2 * p[3] * q[3] * d[3])^2), # 0.9216
    # We don't have HWE, so the above does not hold
    myVarD,
    popVar(gp$dd[, 3, drop = FALSE])[1, 1], # 0.4196571
    # Due to deviation from HWE
    tolerance = 1e-6
  )
  expect_equal(myVarD, unname(gp$varD[3, 3]), tolerance = 1e-6) # 0.4196571

  # ---- Imprinting genetic variance (under Hardy-Weinberg equilibrium) ----

  # No imprinting in loci 1-3
  # expect_equal(gp$varI[1, 1], 0, tolerance = 1e-6)
  # expect_equal(
  #   gp$varI[1, 1],
  #   popVar(gp$i?[, 1, drop = FALSE])[1, 1],
  #   tolerance = 1e-6
  # )
  # expect_equal(gp$varI[2, 2], 0, tolerance = 1e-6)
  # expect_equal(
  #   gp$varI[2, 2],
  #   popVar(gp$i?[, 2, drop = FALSE])[1, 1],
  #   tolerance = 1e-6
  # )
  # expect_equal(gp$varI[3, 3], 0, tolerance = 1e-6)
  # expect_equal(
  #   gp$varI[3, 3],
  #   popVar(gp$i?[, 3, drop = FALSE])[1, 1],
  #   tolerance = 1e-6
  # )

  # ---- Imprinting genetic variance (actual) ----

  # No imprinting in loci 1-3
  # expect_equal(gp$varI[1, 1], 0, tolerance = 1e-6)
  # expect_equal(
  #   gp$varI[1, 1],
  #   popVar(gp$i?[, 1, drop = FALSE])[1, 1],
  #   tolerance = 1e-6
  # )
  # expect_equal(gp$varI[2, 2], 0, tolerance = 1e-6)
  # expect_equal(
  #   gp$varI[2, 2],
  #   popVar(gp$i?[, 2, drop = FALSE])[1, 1],
  #   tolerance = 1e-6
  # )
  # expect_equal(gp$varI[3, 3], 0, tolerance = 1e-6)
  # expect_equal(
  #   gp$varI[3, 3],
  #   popVar(gp$i?[, 3, drop = FALSE])[1, 1],
  #   tolerance = 1e-6
  # )

  # ---- Additive-by-additive genetic variance (under Hardy-Weinberg equilibrium) ----

  # No AxA epistasis in loci 1-3
  expect_equal(gp$varAA[1, 1], 0, tolerance = 1e-6)
  expect_equal(
    gp$varAA[1, 1],
    popVar(gp$aa[, 1, drop = FALSE])[1, 1],
    tolerance = 1e-6
  )
  expect_equal(gp$varAA[2, 2], 0, tolerance = 1e-6)
  expect_equal(
    gp$varAA[2, 2],
    popVar(gp$aa[, 2, drop = FALSE])[1, 1],
    tolerance = 1e-6
  )
  expect_equal(gp$varAA[3, 3], 0, tolerance = 1e-6)
  expect_equal(
    gp$varAA[3, 3],
    popVar(gp$aa[, 3, drop = FALSE])[1, 1],
    tolerance = 1e-6
  )

  # ---- Additive-by-additive genetic variance (actual) ----

  # No AxA epistasis in loci 1-3
  expect_equal(gp$varAA[1, 1], 0, tolerance = 1e-6)
  expect_equal(
    gp$varAA[1, 1],
    popVar(gp$aa[, 1, drop = FALSE])[1, 1],
    tolerance = 1e-6
  )
  expect_equal(gp$varAA[2, 2], 0, tolerance = 1e-6)
  expect_equal(
    gp$varAA[2, 2],
    popVar(gp$aa[, 2, drop = FALSE])[1, 1],
    tolerance = 1e-6
  )
  expect_equal(gp$varAA[3, 3], 0, tolerance = 1e-6)
  expect_equal(
    gp$varAA[3, 3],
    popVar(gp$aa[, 3, drop = FALSE])[1, 1],
    tolerance = 1e-6
  )

  # ---- Additive genic variance (under Hardy-Weinberg) ----

  # TODO: Add genic covariance between traits when there are multiple traits?

  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/136
  # Falconer shows genic variance because he shows one locus only
  theoVarA = unname(2 * p[1] * q[1] * gp$alpha_HW[[1]][1, 1]^2)
  expect_equal(
    theoVarA,
    popVar(gp$bv_HW[, 1, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 1.0368
  expect_equal(theoVarA, unname(gp$genicVarA[1]), tolerance = 1e-6) # 1.0368
  theoVarA = unname(2 * p[2] * q[2] * gp$alpha_HW[[2]][1, 1]^2)
  expect_equal(
    theoVarA,
    popVar(gp$bv_HW[, 2, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 6.2208
  expect_equal(theoVarA, unname(gp$genicVarA[2]), tolerance = 1e-6) # 6.2208
  myVarA = sum(gp$bv_HW[genoLoc3, 3]^2 * c(Q_HW[3], H_HW[3], P_HW[3])) # 6.2208
  expect_equal(
    unname(2 * p[3] * q[3] * gp$alpha_HW[[3]][1, 1]^2), # 6.2208
    myVarA,
    # popVar(gp$bv_HW[, 3, drop = FALSE])[1, 1], # 10.8864
    # Due to deviation from HWE
    tolerance = 1e-6
  )
  expect_equal(myVarA, unname(gp$genicVarA[3]), tolerance = 1e-6) # 6.2208

  # ---- Additive genic variance (actual) ----

  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/136
  # Falconer shows genic variance because he shows one locus only
  theoVarA = unname(2 * p[1] * q[1] * gp$alpha[[1]][1, 1]^2 * (1 + F[1]))
  expect_equal(
    theoVarA,
    popVar(gp$bv[, 1, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 1.0368
  expect_equal(theoVarA, unname(gp$genicVarA[1]), tolerance = 1e-6) # 1.0368
  theoVarA = unname(2 * p[2] * q[2] * gp$alpha[[2]][1, 1]^2 * (1 + F[2]))
  expect_equal(
    theoVarA,
    popVar(gp$bv[, 2, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 6.2208
  expect_equal(theoVarA, unname(gp$genicVarA[2]), tolerance = 1e-6) # 6.2208

  # NOTE: gp$genicVarA is not the actual genic variance, but HWE genic variance;
  #       this deviates from other variances that are actual variances
  #       --> see covA_HW below
  # TODO: Consider changing gp$genicVarA to actual genic variance to match
  #       other variances, for example, VarA is actual genetic variance,
  #       not HWE genetic variance
  myVarA = sum(gp$bv[genoLoc3, 3]^2 * c(Q[3], H[3], P[3])) # 13.05874
  expect_equal(
    # unname(2 * p[3] * q[3] * gp$alpha[[3]][1, 1]^2), # 7.462139
    # unname(2 * p[3] * q[3] * gp$alpha_HW[[3]][1, 1]^2), # 6.2208
    unname(2 * p[3] * q[3] * gp$alpha[[3]][1, 1]^2 * (1 + F[3])), # 13.05874
    myVarA,
    # popVar(gp$bv[, 3, drop = FALSE])[1, 1], # 13.05874
    # Due to deviation from HWE
    tolerance = 1e-6
  )
  # FAILS expect_equal(myVarA, # 13.05874
  #                    unname(gp$genicVarA[3]), tolerance = 1e-6) # 6.2208
  expect_equal(
    unname(2 * p[3] * q[3] * gp$alpha_HW[[3]][1, 1]^2),
    unname(gp$genicVarA[3]),
    tolerance = 1e-6
  ) # 6.2208

  # TODO: why is covA_HW called "additive covariances due to non-random mating"
  #       why "covariance" in particular?
  expect_equal(unname(gp$covA_HW[1]), 0, tolerance = 1e-6) # 0
  expect_equal(unname(gp$covA_HW[2]), 0, tolerance = 1e-6) # 0
  myGenicVarAHWE_vs_actual = unname(
    2 * p[3] * q[3] * gp$alpha[[3]][1, 1]^2 * (1 + F[3])
  ) -
    unname(2 * p[3] * q[3] * gp$alpha_HW[[3]][1, 1]^2)
  # 13.05874 - 6.2208 = 6.83794
  expect_equal(
    myGenicVarAHWE_vs_actual,
    unname(gp$covA_HW[3]),
    tolerance = 1e-6
  ) # 6.83794

  # ---- Dominance genic variance (under Hardy-Weinberg) ----

  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/136
  # Falconer shows genic variance because he shows one locus only
  theoVarD = unname((2 * p[1] * q[1] * d[1])^2)
  expect_equal(
    theoVarD,
    popVar(gp$dd_HW[, 1, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 0.1296
  expect_equal(theoVarD, unname(gp$genicVarD[1]), tolerance = 1e-6) # 0.1296
  theoVarD = unname((2 * p[2] * q[2] * d[2])^2)
  expect_equal(
    theoVarD,
    popVar(gp$dd_HW[, 2, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 0.9216
  expect_equal(theoVarD, unname(gp$genicVarD[2]), tolerance = 1e-6) # 0.9216
  theoVarD = unname((2 * p[3] * q[3] * d[3])^2) # 0.9216
  myVarD = sum(gp$dd_HW[genoLoc3, 3]^2 * c(Q_HW[3], H_HW[3], P_HW[3])) # 0.9216
  expect_equal(
    unname((2 * p[3] * q[3] * d[3])^2), # 0.9216
    myVarD, # 0.9216
    # popVar(gp$dd_HW[, 3, drop = FALSE])[1, 1], # 0.5184
    # Due to deviation from HWE
    tolerance = 1e-6
  )
  expect_equal(myVarD, unname(gp$genicVarD[3]), tolerance = 1e-6) # 0.9216

  # ---- Dominance genic variance (actual) ----

  # Falconer (1996) https://archive.org/details/introductiontoqu0000falc/page/136
  # Falconer shows genic variance because he shows one locus only
  theoVarD = unname((2 * p[1] * q[1] * d[1])^2)
  expect_equal(
    theoVarD,
    popVar(gp$dd[, 1, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 0.1296
  expect_equal(theoVarD, unname(gp$genicVarD[1]), tolerance = 1e-6) # 0.1296
  theoVarD = unname((2 * p[2] * q[2] * d[2])^2)
  expect_equal(
    theoVarD,
    popVar(gp$dd[, 2, drop = FALSE])[1, 1],
    tolerance = 1e-6
  ) # 0.9216
  expect_equal(theoVarD, unname(gp$genicVarD[2]), tolerance = 1e-6) # 0.9216
  theoVarD = unname((2 * p[3] * q[3] * d[3])^2) # 0.9216
  myVarD = sum(gp$dd[genoLoc3, 3]^2 * c(Q[3], H[3], P[3])) # 0.4196571
  expect_equal(
    # theoVarD, # 0.9216
    # We don't have HWE, so the above does not hold
    popVar(gp$dd[, 3, drop = FALSE])[1, 1], # 0.4196571
    myVarD,
    tolerance = 1e-6
  )
  expect_equal(myVarD, unname(gp$genicVarD[3]), tolerance = 1e-6) # 0.4196571
  TODO
  NEXT:model
  the
  dominance
  work
  after
  the
  additive
  work

  # NOTE: gp$genicVarA is not the actual genic variance, but HWE genic variance;
  #       this deviates from other variances that are actual variances
  #       --> see covA_HW below
  # TODO: Consider changing gp$genicVarA to actual genic variance to match
  #       other variances, for example, VarA is actual genetic variance,
  #       not HWE genetic variance
  myVarA = sum(gp$bv[genoLoc3, 3]^2 * c(Q[3], H[3], P[3])) # 13.05874
  expect_equal(
    # unname(2 * p[3] * q[3] * gp$alpha[[3]][1, 1]^2), # 7.462139
    # unname(2 * p[3] * q[3] * gp$alpha_HW[[3]][1, 1]^2), # 6.2208
    unname(2 * p[3] * q[3] * gp$alpha[[3]][1, 1]^2 * (1 + F[3])), # 13.05874
    myVarA,
    # popVar(gp$bv[, 3, drop = FALSE])[1, 1], # 13.05874
    # Due to deviation from HWE
    tolerance = 1e-6
  )
  # FAILS expect_equal(myVarA, # 13.05874
  #                    unname(gp$genicVarA[3]), tolerance = 1e-6) # 6.2208
  expect_equal(
    unname(2 * p[3] * q[3] * gp$alpha_HW[[3]][1, 1]^2),
    unname(gp$genicVarA[3]),
    tolerance = 1e-6
  ) # 6.2208

  # TODO: why is covA_HW called "additive covariances due to non-random mating"
  #       why "covariance" in particular?
  expect_equal(unname(gp$covA_HW[1]), 0, tolerance = 1e-6) # 0
  expect_equal(unname(gp$covA_HW[2]), 0, tolerance = 1e-6) # 0
  myGenicVarAHWE_vs_actual = unname(
    2 * p[3] * q[3] * gp$alpha[[3]][1, 1]^2 * (1 + F[3])
  ) -
    unname(2 * p[3] * q[3] * gp$alpha_HW[[3]][1, 1]^2)
  # 13.05874 - 6.2208 = 6.83794
  expect_equal(
    myGenicVarAHWE_vs_actual,
    unname(gp$covA_HW[3]),
    tolerance = 1e-6
  ) # 6.83794

  # ---- Additive-by-additive genic variance (under Hardy-Weinberg) ----

  # TODO: genicVarAA
  # TODO: covAA_HW

  # ---- Additive-by-additive genic variance (actual) ----

  # TODO: genicVarAA
  # TODO: covAA_HW

  # ---- (Total) Genic variance (under Hardy-Weinberg) ----

  # TODO
  # TODO: genicVarG

  # ---- (Total) Genic variance (actual) ----

  # TODO

  # ---- TODO (under Hardy-Weinberg) ----

  # ---- TODO (actual) ----

  # TODO: covA_L
  # TODO: covD_L
  # TODO: covAA_L
  # TODO: covAA_L
  # TODO: covAD_L
  # TODO: covAAA_L
  # TODO: covDAA_L
  # TODO: covG_L
})
