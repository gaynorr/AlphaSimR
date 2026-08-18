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

test_that("meanG_meanP_meanGPop_meanPPop", {

  founderPop = quickHaplo(nInd=16, nChr=1, segSites=20)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10)
  SP$setVarE(h2=0.5)
  SP$addSnpChip(10)

  pop = newPop(founderPop, simParam=SP)
  ans = RRBLUP(pop, simParam=SP)
  pop = setEBV(pop, ans, simParam=SP)

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

  # Max depth of mp is 3, so level 4 and level 3 return the same output
  mnGV3 = meanGPop(mp, level = 3)
  expect_identical(mnGV3, meanGPop(mp, level = 4))
  mnP3 = meanPPop(mp, level = 3)
  expect_identical(mnP3, meanPPop(mp, level = 4))

  # Level 3: Get meanG and meanP from each independent pop
  calcGV <- meanG(mp, simplify=TRUE, level=1)
  expect_equal(mnGV3, calcGV)
  calcP <- meanP(mp, simplify=TRUE, level=1)
  expect_equal(mnP3, calcP)
  expect_equal(
    meanEBV(mp[1:2], simplify = FALSE),
    lapply(mp@pops[1:2], \(x) lapply(x@pops, \(y) lapply(y@pops, \(z) colMeans(z@ebv))))
  )

  # Level 2: Start calculating mean of means from the previous level
  mnGV2 <- meanGPop(mp, level = 2)
  expect_equal(
    rbind(
      colMeans(meanG(mp$A$e, simplify = TRUE)),
      colMeans(meanG(mp$A$f, simplify = TRUE)),
      colMeans(meanG(mp$B$e, simplify = TRUE)),
      colMeans(meanG(mp$B$f, simplify = TRUE))
    ),
    mnGV2[1:4,,drop=FALSE]
  )
  mnP2 <- meanPPop(mp, level = 2)
  expect_equal(
    rbind(
      colMeans(meanP(mp$A$e, simplify = TRUE)),
      colMeans(meanP(mp$A$f, simplify = TRUE)),
      colMeans(meanP(mp$B$e, simplify = TRUE)),
      colMeans(meanP(mp$B$f, simplify = TRUE))
    ),
    mnP2[1:4,,drop=FALSE]
  )

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
