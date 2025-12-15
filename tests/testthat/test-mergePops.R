# fmt: skip file

context("mergePops")

test_that("mergeMultiPops_and_flattenMultiPop", {

  founderPop = quickHaplo(nInd=100, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10, mean = c(0, 0), var = c(1, 1))
  SP$setVarE(h2=rep(0.5, 2))
  pop = newPop(founderPop, simParam=SP)

  # A pop object is returned unchanged
  expect_identical(pop, mergeMultiPops(pop))
  expect_identical(pop, mergeMultiPops(pop, level = 1))

  # A flat MultiPop
  mp1 = newMultiPop(pop[1:5], pop[6:10])
  # mergePops and mergeMultiPops do the same on MultiPop objects
  expect_identical(mergePops(mp1), mergeMultiPops(mp1))
  # From level 1 of this object, the same input is returned
  expect_identical(mp1, mergeMultiPops(mp1, level = 1))
  expect_identical(mp1@pops[[1]], mergeMultiPops(mp1, level = 1)@pops[[1]])
  expect_identical(mp1@pops[[2]], mergeMultiPops(mp1, level = 1)@pops[[2]])

  # MultiPop with 1 nested object
  mp2 = newMultiPop(pop[1:20], pop[21:40],
                    newMultiPop(pop[41:60],
                                newMultiPop(pop[61:80], pop[81:100])))

  # mergePops doesn't handle nested MultiPop objects
  expect_error(mergePops(mp2), 'all(classes == "Pop") is not TRUE', fixed = TRUE)
  # Use mergePops(list(pop)) instead of just pop because of differences
  #   beteween attributes(pop@ebv) and attributes(mergeMultiPops(mp2)@ebv)
  #   caused by rbind()
  expect_identical(mergePops(list(pop)), mergeMultiPops(mp2))
  expect_identical(pop[1:20],
                   mergeMultiPops(mp2, level = 1)@pops[[1]])
  expect_identical(pop[21:40],
                   mergeMultiPops(mp2, level = 1)@pops[[2]])
  expect_identical(mergePops(list(pop[41:100])),
                   mergeMultiPops(mp2, level = 1)@pops[[3]])
  expect_identical(mp2@pops[[3]][[1]],
                   mergeMultiPops(mp2, level = 2)@pops[[3]][[1]])
  expect_identical(mergePops(list(pop[61:100])),
                   mergeMultiPops(mp2, level = 2)@pops[[3]][[2]])
  expect_identical(mergeMultiPops(mp2, level = 3), mp2)

  # MultiPop with multiple nested objects
  mp3 = newMultiPop(pop[1:20],
                    newMultiPop(pop[21:30],
                                newMultiPop(pop[31:35], pop[36:40])),
                    newMultiPop(pop[41:60],
                                newMultiPop(pop[61:80], pop[81:100])))

  # mergePops doesn't handle nested MultiPop objects
  expect_error(mergePops(mp3), 'all(classes == "Pop") is not TRUE', fixed = TRUE)
  # Use mergePops(list(pop)) instead of just pop because of differences
  #   beteween attributes(pop@ebv) and attributes(mergeMultiPops(mp3)@ebv)
  #   caused by rbind() on empty matrices
  expect_identical(mergePops(list(pop)), mergeMultiPops(mp3))
  expect_identical(pop[1:20],
                   mergeMultiPops(mp3, level = 1)@pops[[1]])
  expect_identical(mergePops(list(pop[21:40])),
                   mergeMultiPops(mp3, level = 1)@pops[[2]])
  expect_identical(mergePops(list(pop[41:100])),
                   mergeMultiPops(mp3, level = 1)@pops[[3]])
  expect_identical(pop[21:30],
                   mergeMultiPops(mp3, level = 2)[[2]][[1]])
  expect_identical(mergePops(list(pop[31:40])),
                   mergeMultiPops(mp3, level = 2)[[2]][[2]])
  expect_identical(pop[41:60],
                   mergeMultiPops(mp3, level = 2)[[3]][[1]])
  expect_identical(mergePops(list(pop[61:100])),
                   mergeMultiPops(mp3, level = 2)[[3]][[2]])
  expect_identical(mergeMultiPops(mp3, level = 3), mp3)

  # Merging multiple objects

  # Merging two pops throws an error
  expect_error(mergeMultiPops(pop[1:15], pop[16:30]),
               "Use mergePops() to merge multiple Pop objects", fixed = TRUE)

  # Merging a pop and a multipop
  expect_identical(mergePops(list(pop[1:10], pop[1:10])),
                   mergeMultiPops(pop[1:10], mp1))
  expect_identical(mergeMultiPops(pop[1:10], mp1, level = 1),
                   c(pop[1:10], mp1))

  # Merging two MultiPops
  expect_identical(mergeMultiPops(mp1, mp2), mergePops(list(pop[1:10], pop)))
  expect_identical(mergeMultiPops(mp1, mp2, level = 1),
                   newMultiPop(pop[1:5], pop[6:10], pop[1:20], pop[21:40],
                               mergePops(list(pop[41:100]))))

  expect_identical(
    mergeMultiPops(mp1, mp2, level = 2),
    newMultiPop(pop[1:5], pop[6:10], pop[1:20], pop[21:40],
                newMultiPop(pop[41:60],
                            mergePops(list(pop[61:80], pop[81:100]))))
  )
  expect_identical(mergeMultiPops(mp1, mp2, level = 3), c(mp1, mp2))

  # Flattening multiPops

  # Nothing to flatten, so the same input should be returned
  expect_identical(flattenMultiPop(mp1, level = 0), mp1)
  expect_identical(flattenMultiPop(mp1, level = 1), mp1)
  expect_identical(flattenMultiPop(mp1, level = 2), mp1)
  tmp = flattenMultiPop(mp1)
  expect_identical(tmp@pops[[1]], pop[1:5])
  expect_identical(tmp@pops[[2]], pop[6:10])
  tmp = flattenMultiPop(mp1, level = 2)
  expect_identical(tmp@pops[[1]], pop[1:5])
  expect_identical(tmp@pops[[2]], pop[6:10])

  # Flatten mp2
  expect_identical(flattenMultiPop(mp2),
                   newMultiPop(pop[1:20], pop[21:40], pop[41:60],
                               pop[61:80], pop[81:100]))
  tmp = flattenMultiPop(mp2, level = 2)
  expect_equal(tmp@pops[[1]], pop[1:20])
  expect_equal(tmp@pops[[2]], pop[21:40])
  expect_equal(tmp@pops[[3]],
               newMultiPop(pop[41:60], pop[61:80], pop[81:100]))
  expect_equal(flattenMultiPop(mp2, level = 3), mp2)

  # Flatten mp3
  expect_equal(flattenMultiPop(mp3),
               newMultiPop(pop[1:20], pop[21:30], pop[31:35], pop[36:40],
                           pop[41:60], pop[61:80], pop[81:100]))
  tmp = flattenMultiPop(mp3, level = 2)
  expect_identical(tmp@pops[[1]], pop[1:20])
  expect_identical(tmp@pops[[2]],
                   newMultiPop(pop[21:30], pop[31:35], pop[36:40]))
  expect_identical(tmp@pops[[3]],
                   newMultiPop(pop[41:60], pop[61:80], pop[81:100]))
  expect_equal(flattenMultiPop(mp3, level = 3), mp3)
})
