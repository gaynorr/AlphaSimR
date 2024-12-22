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

  expect_error(parentAverage(pop2, mothers = pop, fathers = pop, use = "x"))
  expect_error(parentAverage(pop2, mothers = pop, fathers = pop, use = "aa"))
  expect_error(parentAverage(pop2, mothers = pop, fathers = pop, use = "dd"))
  expect_error(parentAverage(pop2, mothers = pop, fathers = pop, use = "id"))

  pa_gv = parentAverage(pop2, mothers = pop, fathers = pop, use = "gv", simParam=SP)
  ms_gv = mendelianSampling(pop2, mothers = pop, fathers = pop, use = "gv", simParam=SP)

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

