# fmt: skip file
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
