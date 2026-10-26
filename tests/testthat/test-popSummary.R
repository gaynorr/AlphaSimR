context("popSummary")

test_that("meanP, meanEBV, varP, and varEBV work", {
  set.seed(826)
  founderPop = quickHaplo(nInd=3, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10, mean=c(0, 0), var=c(1, 1))
  pop = newPop(founderPop, simParam=SP)
  pop@pheno = cbind(first=c(2, 4, 9), second=c(8, 2, 5))
  pop@ebv = cbind(first=c(1, 3, 5), second=c(2, 6, 4))

  expect_equal(meanP(pop), c(first=5, second=5))
  expect_equal(meanEBV(pop), c(first=3, second=4))

  # The off-diagonal entries and divisor distinguish covariance from sample variance.
  expectedVarP = matrix(c(26, -6, -6, 18) / 3, nrow=2,
                    dimnames=list(c("first", "second"), c("first", "second")))
  expectedVarEBV = matrix(c(8, 4, 4, 8) / 3, nrow=2,
                    dimnames=list(c("first", "second"), c("first", "second")))
  expect_equal(varP(pop), expectedVarP)
  expect_equal(varEBV(pop), expectedVarEBV)
  expect_equal(meanP(pop[1]), c(first=2, second=8))
  expect_equal(meanEBV(pop[1]), c(first=1, second=2))
  expect_equal(varP(pop[1]), expectedVarP * 0)
  expect_equal(varEBV(pop[1]), expectedVarEBV * 0)

  SP1 = SimParam$new(founderPop)
  SP1$nThreads = 1L
  SP1$addTraitA(10)
  oneTrait = newPop(founderPop, simParam=SP1)
  oneTrait@pheno = pop@pheno[, "first", drop=FALSE]
  oneTrait@ebv = pop@ebv[, "first", drop=FALSE]
  expect_equal(meanP(oneTrait), c(first=5))
  expect_equal(meanEBV(oneTrait), c(first=3))
  expect_equal(varP(oneTrait), expectedVarP["first", "first", drop=FALSE])
  expect_equal(varEBV(oneTrait), expectedVarEBV["first", "first", drop=FALSE])
  expect_equal(meanP(oneTrait[1]), c(first=2))
  expect_equal(meanEBV(oneTrait[1]), c(first=1))
  expect_equal(varP(oneTrait[1]), expectedVarP["first", "first", drop=FALSE] * 0)
  expect_equal(varEBV(oneTrait[1]), expectedVarEBV["first", "first", drop=FALSE] * 0)
})

test_that("parentAverage and mendelianSampling work",{
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
  # Different EBVs expose accidental use of the genetic-value slot.
  pop@ebv = pop@gv + 10
  pop2@ebv = pop2@gv + 20

  for (fun in list(parentAverage, mendelianSampling)) {
    for (args in list(list(), list(mothers=pop), list(fathers=pop))) {
      expect_error(do.call(fun, c(list(pop=pop2, simParam=SP), args)),
                   "must provide either 'parents' or both 'mothers' and 'fathers'!",
                   fixed=TRUE)
    }
    expect_error(fun(pop2, parents=pop[1:2], simParam=SP),
                 "some parents/mothers not found!", fixed=TRUE)
    expect_error(fun(pop2, parents=pop[c(1, 3)], simParam=SP),
                 "some parents/fathers not found!", fixed=TRUE)
    expect_error(fun(pop2, mothers=pop[1], fathers=pop[2], simParam=SP),
                 "some parents/mothers not found!", fixed=TRUE)
    expect_error(fun(pop2, mothers=pop[c(1, 3)], fathers=pop[1], simParam=SP),
                 "some parents/fathers not found!", fixed=TRUE)
    for (use in c("x", "aa", "dd", "id")) {
      expect_error(fun(pop2, mothers=pop, fathers=pop, use=use, simParam=SP),
                   "use must be one of 'gv', 'ebv', or 'pheno'!", fixed=TRUE)
    }
  }

  # Parent matching must use IDs, regardless of row order or separate parent subsets
  pa_gv = parentAverage(pop2, parents = pop[c(3, 1, 2)], use = "gv", simParam=SP)
  ms_gv = mendelianSampling(pop2, parents = pop[c(3, 1, 2)], use = "gv", simParam=SP)
  pa_gv2 = parentAverage(pop2, mothers = pop[c(3, 1)], fathers = pop[2], use = "gv", simParam=SP)
  ms_gv2 = mendelianSampling(pop2, mothers = pop[c(3, 1)], fathers = pop[2], use = "gv", simParam=SP)

  expect_equal(pa_gv, pa_gv2)
  expect_equal(ms_gv, ms_gv2)

  # One offspring must still return a matrix, for either way of supplying parents
  expect_equal(parentAverage(pop2[1], parents=pop, simParam=SP), pa_gv[1, , drop=FALSE])
  expect_equal(mendelianSampling(pop2[1], parents=pop, simParam=SP), ms_gv[1, , drop=FALSE])
  expect_equal(parentAverage(pop2[1], mothers=pop, fathers=pop, simParam=SP), pa_gv[1, , drop=FALSE])
  expect_equal(mendelianSampling(pop2[1], mothers=pop, fathers=pop, simParam=SP), ms_gv[1, , drop=FALSE])

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

  for (fun in list(parentAverage, mendelianSampling)) {
    expect_error(fun(pop2, parents = pop, use = "bv", simParam=SP),
                 "use must be one of 'gv', 'ebv', or 'pheno'!", fixed=TRUE)
    expect_error(fun(pop2, mothers = pop, fathers = pop, use = "bv", simParam=SP),
                 "use must be one of 'gv', 'ebv', or 'pheno'!", fixed=TRUE)
  }

  expect_identical(parentAverage(pop2, parents = pop, use = "gv", simParam=SP, nThreads = 2), pa_gv)
  expect_identical(mendelianSampling(pop2, parents = pop, use = "gv", simParam=SP, nThreads = 2), ms_gv)

  local({
    hadSP = exists("SP", envir=.GlobalEnv, inherits=FALSE)
    if (hadSP) oldSP = get("SP", envir=.GlobalEnv)
    on.exit({
      if (hadSP) assign("SP", oldSP, envir=.GlobalEnv)
      else rm("SP", envir=.GlobalEnv)
    })
    # Exercise the global default without retaining changes to the user's SP
    assign("SP", SP, envir=.GlobalEnv)
    expect_identical(genParam(pop2), genParam(pop2, simParam=SP))
    expect_identical(parentAverage(pop2, parents=pop),
                     parentAverage(pop2, parents=pop, simParam=SP))
    expect_identical(mendelianSampling(pop2, parents=pop),
                     mendelianSampling(pop2, parents=pop, simParam=SP))
  })

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
      function(x) sample(rep(letters[5:6], length.out = length(x))),
      function(x) sample(rep(letters[7:8], length.out = length(x)))
    )
  )
  idx = sample(1:16, 8)
  mp$C = splitPop(pop[idx], function(x) sample(rep(letters[5:6], length.out = length(x))))
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
    lapply(mp@pops[1:2], function(x) lapply(x@pops, function(y) lapply(y@pops, function(z) colMeans(z@ebv))))
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

  # Duplicate MultiPop names preserve source rows
    mp = newMultiPop(
      dup = newMultiPop(
        left = pop[1:2],
        right = pop[3:4]
      ),
      dup = newMultiPop(
        left = pop[5:6],
        right = pop[7:8]
      )
    )

    result = meanGPop(mp, level = 1)
    source = attr(result, "source")

    expect_equal(nrow(result), 2L)
    expect_equal(nrow(source), 2L)
    expect_false(identical(source[1, ], source[2, ]))
  
  # A MultiPop from independently created pop objects
  founderPop2 = quickHaplo(nInd = 4, nChr = 1, segSites = 10)
  SP2 = SimParam$new(founderPop2)
  SP2$nThreads = 1L
  SP2$addTraitA(10, mean = 10, var = 1, name = 'A1')
  SP2$setVarE(h2 = 0.5)

  pop2 = newPop(founderPop2, simParam = SP2)

  chimera_mp = newMultiPop(
    dup = newMultiPop(
      left = pop[1:2],
      right = pop[3:4]
    ),
    dup = newMultiPop(
      left = pop2[1:2],
      right = pop2[3:4]
    )
  )

  # Duplicate MultiPop names preserve source rows
    mp = newMultiPop(
      dup = newMultiPop(
        left = pop[1:2],
        right = pop[3:4]
      ),
      dup = newMultiPop(
        left = pop[5:6],
        right = pop[7:8]
      )
    )

    result = meanGPop(mp, level = 1)
    source = attr(result, "source")

    expect_equal(nrow(result), 2L)
    expect_equal(nrow(source), 2L)
    expect_false(identical(source[1, ], source[2, ]))
  
  # A MultiPop from independently created pop objects
  founderPop2 = quickHaplo(nInd = 4, nChr = 1, segSites = 10)
  SP2 = SimParam$new(founderPop2)
  SP2$nThreads = 1L
  SP2$addTraitA(10, mean = 10, var = 1, name = 'A1')
  SP2$setVarE(h2 = 0.5)

  pop2 = newPop(founderPop2, simParam = SP2)

  chimera_mp = newMultiPop(
    dup = newMultiPop(
      left = pop[1:2],
      right = pop[3:4]
    ),
    dup = newMultiPop(
      left = pop2[1:2],
      right = pop2[3:4]
    )
  )

  # Expected errors and warnings
  expect_warning(
    gv(chimera_mp, simplify = TRUE),
    "Some values returned by FUN do not have consistent column names. Returning list output."
  )
  expect_warning(
    gv(chimera_mp, simplify = TRUE),
    "Some values returned by FUN do not have consistent column names. Returning list output."
  )
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
    meanGPop(mp, level = 0.5),
    "`level` must be a single non-NA integer value."
  )
  expect_error(
    meanGPop(mp, level = 1:2),
    "`level` must be a single non-NA integer value."
  )
  expect_error(
    meanPPop(mp, level = 0.5),
    "`level` must be a single non-NA integer value."
  )
  expect_error(
    meanPPop(mp, level = 1:2),
    "`level` must be a single non-NA integer value."
  )
  expect_error(
    meanGPop(chimera_mp),
    "Population values do not have consistent trait names."
  )
  expect_error(
    meanPPop(chimera_mp),
    "Population values do not have consistent trait names."
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
})

test_that("HYBRIDPOP summaries accept a HybridPop", {
  founderPop = quickHaplo(nInd=6, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10, mean=c(0, 0), var=c(1, 1))
  SP$setVarE(h2=c(0.6, 0.4))
  pop = newPop(founderPop, simParam=SP)
  hybrid = hybridCross(pop[1:3], pop[4:6], returnHybridPop=TRUE,
                       simParam=SP)
  expect_true(isHybridPop(hybrid))

  expect_equal(meanG(hybrid), colMeans(hybrid@gv))
  expect_equal(meanP(hybrid), colMeans(hybrid@pheno))
  expect_equal(gv(hybrid), hybrid@gv)
  expect_equal(pheno(hybrid), hybrid@pheno)
  expect_equivalent(varG(hybrid), popVar(hybrid@gv))
  expect_equivalent(varP(hybrid), popVar(hybrid@pheno))
})
