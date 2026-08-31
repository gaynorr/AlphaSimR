# fmt: skip file
context("popSummary")

test_that("parentAverage_and_mendelianSampling_work",{
  founderPop = quickHaplo(nInd=3, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$addTraitAD(10, mean=c(0, 0), var=c(1, 1), meanDD=c(0, 0.5))
  SP$setVarE(varE=c(0.5, 0.5))
  SP$nThreads = 1L
  pop = newPop(founderPop, simParam=SP)
  pop2 = makeCross(pop, crossPlan = matrix(data=c(1, 2,
                                                  3, 2),
                                           byrow=TRUE, ncol=2),
                   nProgeny=2, simParam=SP)
  pop@ebv = pop@gv
  pop2@ebv = pop2@gv

  expect_error(parentAverage(pop2, mothers = pop))
  expect_error(parentAverage(pop2, fathers = pop))
  expect_error(parentAverage(pop2, parents = pop[1:2]))
  expect_error(parentAverage(pop2, mothers = pop, fathers = pop, use = "x"))
  expect_error(parentAverage(pop2, mothers = pop, fathers = pop, use = "aa"))
  expect_error(parentAverage(pop2, mothers = pop, fathers = pop, use = "dd"))
  expect_error(parentAverage(pop2, mothers = pop, fathers = pop, use = "id"))

  pa_gv = parentAverage(pop2, parents = pop, use = "gv", simParam=SP)
  ms_gv = mendelianSampling(pop2, parents = pop, use = "gv", simParam=SP)
  pa_gv2 = parentAverage(pop2, mothers = pop, fathers = pop, use = "gv", simParam=SP)
  ms_gv2 = mendelianSampling(pop2, mothers = pop, fathers = pop, use = "gv", simParam=SP)

  expect_equal(pa_gv, pa_gv2)
  expect_equal(ms_gv, ms_gv2)

  expect_equal(pa_gv[1, ], pa_gv[2, ])
  expect_equal(pa_gv[3, ], pa_gv[4, ])

  expect_equal(pa_gv[1, ], 0.5 * (pop@gv[1, ] + pop@gv[2, ]))
  expect_equal(pa_gv[2, ], 0.5 * (pop@gv[1, ] + pop@gv[2, ]))
  expect_equal(pa_gv[3, ], 0.5 * (pop@gv[3, ] + pop@gv[2, ]))
  expect_equal(pa_gv[4, ], 0.5 * (pop@gv[3, ] + pop@gv[2, ]))

  expect_equal(ms_gv[1, ], pop2@gv[1, ] - 0.5 * (pop@gv[1, ] + pop@gv[2, ]))
  expect_equal(ms_gv[2, ], pop2@gv[2, ] - 0.5 * (pop@gv[1, ] + pop@gv[2, ]))
  expect_equal(ms_gv[3, ], pop2@gv[3, ] - 0.5 * (pop@gv[3, ] + pop@gv[2, ]))
  expect_equal(ms_gv[4, ], pop2@gv[4, ] - 0.5 * (pop@gv[3, ] + pop@gv[2, ]))

  pa_ebv = parentAverage(pop2, mothers = pop, fathers = pop, use = "ebv", simParam=SP)
  ms_ebv = mendelianSampling(pop2, mothers = pop, fathers = pop, use = "ebv", simParam=SP)
  expect_equal(pa_ebv[1, ], pa_ebv[2, ])
  expect_equal(ms_ebv[1, ], pop2@ebv[1, ] - 0.5 * (pop@ebv[1, ] + pop@ebv[2, ]))

  pa_pheno = parentAverage(pop2, mothers = pop, fathers = pop, use = "pheno", simParam=SP)
  ms_pheno = mendelianSampling(pop2, mothers = pop, fathers = pop, use = "pheno", simParam=SP)
  expect_equal(pa_pheno[1, ], pa_pheno[2, ])
  expect_equal(ms_pheno[1, ], pop2@pheno[1, ] - 0.5 * (pop@pheno[1, ] + pop@pheno[2, ]))

  pa_bv = parentAverage(pop2, mothers = pop, fathers = pop, use = "bv", simParam=SP)
  ms_bv = mendelianSampling(pop2, mothers = pop, fathers = pop, use = "bv", simParam=SP)
  expect_equal(pa_bv[1, ], pa_bv[2, ])
  expect_equal(ms_bv[1, ], bv(pop2, simParam=SP)[1, ] - 0.5 * (bv(pop, simParam=SP)[1, ] + bv(pop, simParam=SP)[2, ]))
})

test_that("slots_means_and_variances", {
  
  founderPop = quickHaplo(nInd=16, nChr=1, segSites=20)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10, mean=c(0,0), var=c(1,1))
  SP$setVarE(h2=c(0.5, 0.2))
  SP$addSnpChip(10)

  pop = newPop(founderPop, simParam=SP)
  expect_equal(nPop(pop), 1L)
  ans = RRBLUP(pop, simParam=SP)
  pop = setEBV(pop, ans, simParam=SP)

  # Check @gv @pheno and @ebv slots
  expect_equal(gv(pop), pop@gv)
  expect_equal(pheno(pop), pop@pheno)
  expect_equal(ebv(pop), pop@ebv)

  # Check means and variances of a single pop
  expect_equal(meanG(pop), colMeans(pop@gv))
  expect_equal(meanP(pop), colMeans(pop@pheno))
  expect_equal(meanEBV(pop), colMeans(pop@ebv))
  expect_equivalent(varG(pop), popVar(pop@gv))
  expect_equivalent(varP(pop), popVar(pop@pheno))
  expect_equivalent(varEBV(pop), popVar(pop@ebv))

  # Create a multipop with complex nesting
  mp = splitPop(
    pop,
    by = list(
      sample(rep(LETTERS[1:2], length.out = pop@nInd)),
      \(x) sample(rep(letters[5:6], length.out = length(x))),
      \(x) sample(rep(letters[7:8], length.out = length(x)))
    )
  )
  idx = sample(1:16, 8)
  mp$C = splitPop(pop[idx], \(x) sample(rep(letters[5:6], length.out = length(x))))
  mp$D = pop[-idx]

  # Verify the number of terminal Pop objects
  expect_equal(nPop(mp), 11L)

  # Check within-Pop variances in a MultiPop
  vrG = varG(mp)
  expect_equal(vrG$A$e$g, varG(mp$A$e$g))
  expect_equal(vrG$A$e$h, varG(mp$A$e$h))
  expect_equal(vrG$A$f$g, varG(mp$A$f$g))
  expect_equal(vrG$A$f$h, varG(mp$A$f$h))
  expect_equal(vrG$B$e$g, varG(mp$B$e$g))
  expect_equal(vrG$B$e$h, varG(mp$B$e$h))
  expect_equal(vrG$B$f$g, varG(mp$B$f$g))
  expect_equal(vrG$B$f$h, varG(mp$B$f$h))
  expect_equal(vrG$C$e, varG(mp$C$e))
  expect_equal(vrG$C$f, varG(mp$C$f))
  expect_equal(vrG$D, varG(mp$D))

  vrP = varP(mp)
  expect_equal(vrP$A$e$g, varP(mp$A$e$g))
  expect_equal(vrP$A$e$h, varP(mp$A$e$h))
  expect_equal(vrP$A$f$g, varP(mp$A$f$g))
  expect_equal(vrP$A$f$h, varP(mp$A$f$h))
  expect_equal(vrP$B$e$g, varP(mp$B$e$g))
  expect_equal(vrP$B$e$h, varP(mp$B$e$h))
  expect_equal(vrP$B$f$g, varP(mp$B$f$g))
  expect_equal(vrP$B$f$h, varP(mp$B$f$h))
  expect_equal(vrP$C$e, varP(mp$C$e))
  expect_equal(vrP$C$f, varP(mp$C$f))
  expect_equal(vrP$D, varP(mp$D))

  vrEBV = varEBV(mp)
  expect_equal(vrEBV$A$e$g, varEBV(mp$A$e$g))
  expect_equal(vrEBV$A$e$h, varEBV(mp$A$e$h))
  expect_equal(vrEBV$A$f$g, varEBV(mp$A$f$g))
  expect_equal(vrEBV$A$f$h, varEBV(mp$A$f$h))
  expect_equal(vrEBV$B$e$g, varEBV(mp$B$e$g))
  expect_equal(vrEBV$B$e$h, varEBV(mp$B$e$h))
  expect_equal(vrEBV$B$f$g, varEBV(mp$B$f$g))
  expect_equal(vrEBV$B$f$h, varEBV(mp$B$f$h))
  expect_equal(vrEBV$C$e, varEBV(mp$C$e))
  expect_equal(vrEBV$C$f, varEBV(mp$C$f))
  expect_equal(vrEBV$D, varEBV(mp$D))

  # Max depth of mp is 3, so level 4 and level 3 return the same output
  mnGV3 = meanGPop(mp, level = 3)
  expect_identical(mnGV3, meanGPop(mp, level = 4))
  mnP3 = meanPPop(mp, level = 3)
  expect_identical(mnP3, meanPPop(mp, level = 4))

  # Level 3: Get meanG and meanP from each independent pop
  calcGV3 = meanG(mp, simplify=TRUE, level=1)
  expect_equal(mnGV3, calcGV3)
  calcP3 = meanP(mp, simplify=TRUE, level=1)
  expect_equal(mnP3, calcP3)
  expect_equal(
    meanEBV(mp[1:2], simplify = FALSE),
    lapply(mp@pops[1:2], \(x) lapply(x@pops, \(y) lapply(y@pops, \(z) colMeans(z@ebv))))
  )

  # Level 3: Calculate between-pop variance
  expect_equivalent(varGPop(mp, level = 3), popVar(mnGV3))
  expect_equivalent(varPPop(mp, level = 3), popVar(mnP3))

  # Level 2: Start calculating mean of means from the previous level
  mnGV2 = meanGPop(mp, level = 2)
  expect_equal(
    rbind(
      colMeans(meanG(mp$A$e, simplify = TRUE)),
      colMeans(meanG(mp$A$f, simplify = TRUE)),
      colMeans(meanG(mp$B$e, simplify = TRUE)),
      colMeans(meanG(mp$B$f, simplify = TRUE))
    ),
    mnGV2[1:4,,drop=FALSE]
  )
  mnP2 = meanPPop(mp, level = 2)
  expect_equal(
    rbind(
      colMeans(meanP(mp$A$e, simplify = TRUE)),
      colMeans(meanP(mp$A$f, simplify = TRUE)),
      colMeans(meanP(mp$B$e, simplify = TRUE)),
      colMeans(meanP(mp$B$f, simplify = TRUE))
    ),
    mnP2[1:4,,drop=FALSE]
  )

  # Level 2: Calculate between-pop variance
  expect_equivalent(varGPop(mp, level = 2), popVar(mnGV2))
  expect_equivalent(varPPop(mp, level = 2), popVar(mnP2))

  # Level 1
  mnGV1 = meanGPop(mp, level = 1)
  expect_equal(
    rbind(
      colMeans(mnGV2[1:2,,drop=FALSE]),
      colMeans(mnGV2[3:4,,drop=FALSE ]),
      colMeans(mnGV2[5:6,,drop=FALSE ]),
      colMeans(mnGV2[7  ,,drop=FALSE])
    ),
    mnGV1[,,drop=FALSE]
  )
  mnP1 = meanPPop(mp, level = 1)
  expect_equal(
    rbind(
      colMeans(mnP2[1:2,,drop=FALSE]),
      colMeans(mnP2[3:4,,drop=FALSE ]),
      colMeans(mnP2[5:6,,drop=FALSE ]),
      colMeans(mnP2[7  ,,drop=FALSE])
    ),
    mnP1[,,drop=FALSE]
  )

  # Level 1: Calculate between-pop variance
  expect_equivalent(varGPop(mp, level = 1), popVar(mnGV1))
  expect_equivalent(varPPop(mp, level = 1), popVar(mnP1))

  # Level 0
  expect_equal(
    meanGPop(mp, level = 0),
    colMeans(mnGV1)
  )
  expect_equal(
    meanPPop(mp, level = 0),
    colMeans(mnP1)
  )

  # Expected errors and warnings
  expect_warning(
    gv(mp, simplify = TRUE, level = 0),
    "`level` should be >= 1. Setting default `level=1`"
  )
  expect_warning(
    pheno(mp, simplify = TRUE, level = 0),
    "`level` should be >= 1. Setting default `level=1`"
  )
  expect_warning(
    ebv(mp, simplify = TRUE, level = 0),
    "`level` should be >= 1. Setting default `level=1`"
  )
  expect_warning(
    meanG(mp, simplify = TRUE, level = 0),
    "`level` should be >= 1. Setting default `level=1`"
  )
  expect_warning(
    meanP(mp, simplify = TRUE, level = 0),
    "`level` should be >= 1. Setting default `level=1`"
  )
  expect_warning(
    meanEBV(mp, simplify = TRUE, level = 0),
    "`level` should be >= 1. Setting default `level=1`"
  )
  expect_warning(
    varGPop(mp, level = 0),
    "Returning a variance-covariance matrix of zeros. You may want to increase the value of `level`."
  )
  expect_warning(
    varPPop(mp, level = 0),
    "Returning a variance-covariance matrix of zeros. You may want to increase the value of `level`."
  )
  expect_error(
    meanGPop(newEmptyPop(ploidy = 2L, simParam = SP)),
    "One of the populations in `x` is empty",
    fixed = TRUE
  )
  expect_error(
    meanPPop(newEmptyPop(ploidy = 2L, simParam = SP)),
    "One of the populations in `x` is empty",
    fixed = TRUE
  )
  expect_error(
    meanGPop(newEmptyMultiPop()),
    "`x` contains no populations.",
    fixed = TRUE
  )
  expect_error(
    meanPPop(newEmptyMultiPop()),
    "`x` contains no populations.",
    fixed = TRUE
  )
  expect_error(
    meanGPop(pop, level = c(NA, NA)),
    "`level` must be a single non-NA integer value.",
    fixed = TRUE
  )
  expect_error(
    meanPPop(pop, level = c(NA, NA)),
    "`level` must be a single non-NA integer value.",
    fixed = TRUE
  )
})