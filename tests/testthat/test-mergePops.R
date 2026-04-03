# fmt: skip file

context("mergePops")

test_that("cPop_and_mergePops", {
  founderPop = quickHaplo(nInd=3, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  pop = newPop(founderPop, simParam=SP)
  expect_identical(c(pop[1:2], pop[3]), pop[1:3])
  expect_identical(mergePops(list(pop[1:2], pop[3])), pop[1:3])

  founderPop = quickHaplo(nInd=3, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10)
  SP$setVarE(h2=0.5)
  expect_identical(c(pop[1:2], pop[3]), pop[1:3])
  expect_identical(mergePops(list(pop[1:2], pop[3])), pop[1:3])

  pop@ebv = pop@pheno
  expect_identical(c(pop[1:2], pop[3]), pop[1:3])
  expect_identical(mergePops(list(pop[1:2], pop[3])), pop[1:3])
})

test_that("cMultiPop_mergeMultiPops_and_flattenMultiPop", {

  founderPop = quickHaplo(nInd=12, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10)
  SP$setVarE(h2=0.5)
  pop = newPop(founderPop, simParam=SP)

  # A single Pop object is returned unchanged
  expect_identical(pop, mergeMultiPops(pop))
  expect_identical(pop, mergeMultiPops(pop, level = 1))
  
  # A flat MultiPop
  mp1 = newMultiPop(pop[1:2], pop[3:5])
  expect_identical(c(mp1), mp1)
  # mergePops and mergeMultiPops do the same on MultiPop objects
  expect_identical(mergePops(mp1), mergeMultiPops(mp1))
  # From level 1 of this object, the same input is returned
  expect_identical(mp1, mergeMultiPops(mp1, level = 1))
  expect_identical(mp1[[1]], mergeMultiPops(mp1, level = 1)[[1]])
  expect_identical(mp1[[2]], mergeMultiPops(mp1, level = 1)[[2]])
  expect_identical(mp1, mergeMultiPops(mp1, level = 2))

  # MultiPop with one nested object
  mp2 = newMultiPop(pop[1:2], pop[3:5],
                    newMultiPop(pop[6:7],
                                newMultiPop(pop[8], pop[9:10])))

  # mergePops doesn't handle nested MultiPop objects
  expect_error(mergePops(mp2), 'all(classes == "Pop") is not TRUE', fixed = TRUE)
  expect_identical(pop[1:10], mergeMultiPops(mp2))
  expect_identical(pop[1:2],
                   mergeMultiPops(mp2, level = 1)[[1]])
  expect_identical(pop[3:5],
                   mergeMultiPops(mp2, level = 1)[[2]])
  expect_identical(mergePops(list(pop[6:10])),
                   mergeMultiPops(mp2, level = 1)[[3]])
  expect_identical(mp2[[3]][[1]],
                   mergeMultiPops(mp2, level = 2)[[3]][[1]])
  expect_identical(mergePops(list(pop[8:10])),
                   mergeMultiPops(mp2, level = 2)[[3]][[2]])
  expect_identical(mergeMultiPops(mp2, level = 3), mp2)

  # MultiPop with multiple nested objects
  mp3 = newMultiPop(pop[1:2],
                    newMultiPop(pop[3:4],
                                newMultiPop(pop[5], pop[6:7])),
                    newMultiPop(pop[8:9],
                                newMultiPop(pop[10:11], pop[12])))

  # mergePops doesn't handle nested MultiPop objects
  expect_error(mergePops(mp3), 'all(classes == "Pop") is not TRUE', fixed = TRUE)
  expect_identical(pop, mergeMultiPops(mp3))
  expect_identical(pop[1:2],
                   mergeMultiPops(mp3, level = 1)[[1]])
  expect_identical(mergePops(list(pop[3:7])),
                   mergeMultiPops(mp3, level = 1)[[2]])
  expect_identical(mergePops(list(pop[8:12])),
                   mergeMultiPops(mp3, level = 1)[[3]])
  expect_identical(pop[3:4],
                   mergeMultiPops(mp3, level = 2)[[2]][[1]])
  expect_identical(mergePops(list(pop[5:7])),
                   mergeMultiPops(mp3, level = 2)[[2]][[2]])
  expect_identical(pop[8:9],
                   mergeMultiPops(mp3, level = 2)[[3]][[1]])
  expect_identical(mergePops(list(pop[10:12])),
                   mergeMultiPops(mp3, level = 2)[[3]][[2]])
  expect_identical(mergeMultiPops(mp3, level = 3), mp3)

  # Merging multiple objects

  # Merging two pops throws an error
  expect_error(mergeMultiPops(pop[1:2], pop[3:4]),
               "Use mergePops() to merge multiple Pop objects", fixed = TRUE)
  # Merging invalid objects throws an error
  expect_error(mergeMultiPops(pop, SP),
               'all(classes == "Pop") is not TRUE', fixed = TRUE)
  expect_error(mergeMultiPops(mp1, SP),
               'is(y, "MultiPop") is not TRUE', fixed = TRUE)
  expect_error(mergeMultiPops(pop, SP, mp1),
               'is(y, "MultiPop") is not TRUE', fixed = TRUE)
  expect_error(mergeMultiPops(mp1, SP, pop),
               'is(y, "MultiPop") is not TRUE', fixed = TRUE)

  # Combining a Pop and a MultiPop gives a MultiPop
  expect_identical(c(mp1, pop[6:10]),
                   newMultiPop(pop[1:2], pop[3:5], pop[6:10]))
  # ... testing the other order to check c() works correctly
  expect_identical(c(pop[6:10], mp1),
                   newMultiPop(pop[6:10], pop[1:2], pop[3:5]))
  expect_identical(c(mp1, pop[6:10], mp2),
                   mergeMultiPops(mp1, pop[6:10], mp2, level = 3))

  # Merging a pop and a multipop
  expect_identical(pop[1:10],
                   mergeMultiPops(mp1, pop[6:10]))
  expect_identical(mergeMultiPops(mp1, pop[6:10], level = 1),
                   c(mp1, pop[6:10]))
  # Note that c(mp1, pop) gives a MultiPop, but doesn't have level control
  expect_identical(c(pop[6:10], mp1),
                   mergeMultiPops(pop[6:10], mp1, level = 1))


  # Merging two MultiPops
  expect_identical(mergeMultiPops(mp1, mp2), c(pop[1:5], pop[1:10]))
  expect_identical(mergeMultiPops(mp1, mp2, level = 1),
                   newMultiPop(pop[1:2], pop[3:5], pop[1:2], pop[3:5], pop[6:10]))

  expect_identical(
    mergeMultiPops(mp1, mp2, level = 2),
    newMultiPop(pop[1:2], pop[3:5], pop[1:2], pop[3:5],
                newMultiPop(pop[6:7], pop[8:10]))
  )
  expect_identical(mergeMultiPops(mp1, mp2, level = 3), c(mp1, mp2))

  # Flattening multiPops

  # Nothing to flatten, so the same input should be returned
  expect_identical(flattenMultiPop(mp1, level = 0), mp1)
  expect_identical(flattenMultiPop(mp1, level = 1), mp1)
  expect_identical(flattenMultiPop(mp1, level = 2), mp1)
  tmp = flattenMultiPop(mp1)
  expect_identical(tmp[[1]], pop[1:2])
  expect_identical(tmp[[2]], pop[3:5])
  tmp = flattenMultiPop(mp1, level = 2)
  expect_identical(tmp[[1]], pop[1:2])
  expect_identical(tmp[[2]], pop[3:5])

  # Flatten mp2
  expect_identical(flattenMultiPop(mp2),
                   newMultiPop(pop[1:2], pop[3:5], pop[6:7],
                               pop[8], pop[9:10]))
  tmp = flattenMultiPop(mp2, level = 2)
  expect_equal(tmp[[1]], pop[1:2])
  expect_equal(tmp[[2]], pop[3:5])
  expect_equal(tmp[[3]],
               newMultiPop(pop[6:7], pop[8], pop[9:10]))
  expect_equal(flattenMultiPop(mp2, level = 3), mp2)

  # Flatten mp3
  expect_equal(flattenMultiPop(mp3),
               newMultiPop(pop[1:2], pop[3:4], pop[5], pop[6:7],
                           pop[8:9], pop[10:11], pop[12]))
  tmp = flattenMultiPop(mp3, level = 2)
  expect_identical(tmp[[1]], pop[1:2])
  expect_identical(tmp[[2]],
                   newMultiPop(pop[3:4], pop[5], pop[6:7]))
  expect_identical(tmp[[3]],
                   newMultiPop(pop[8:9], pop[10:11], pop[12]))
  expect_equal(flattenMultiPop(mp3, level = 3), mp3)
})
