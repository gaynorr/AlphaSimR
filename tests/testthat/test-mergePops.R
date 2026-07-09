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
  expect_identical(mergePops(list(pop[1:2], pop[3], NULL)), pop[1:3])
  expect_identical(mergePops(newMultiPop(pop[1:2], newMultiPop(pop[3]))), pop[1:3])
  
  pop@ebv = pop@pheno
  expect_identical(c(pop[1:2], pop[3]), pop[1:3])
  expect_identical(mergePops(list(pop[1:2], pop[3])), pop[1:3])
  
  pop3 = pop[3]
  pop3@ebv = matrix(1:2, ncol=2, dimnames=list(NULL, c("Trait1", "Trait2")))
  expect_warning(c(pop[1:2], pop3),
                 "Populations have different numbers of EBV columns; EBVs removed!", fixed = TRUE)
  expect_warning(mergePops(list(pop[1:2], pop3)),
                 "Populations have different numbers of EBV columns; EBVs removed!", fixed = TRUE)
})

test_that("cMultiPop_mergeMultiPops_and_flattenMultiPop", {

  founderPop = quickHaplo(nInd=12, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10)
  SP$setVarE(h2=0.5)
  
  # Non-Pop and non-MultiPop objects throw an error
  expect_error(mergeMultiPops(1:5), 
               "One or more objects are not of Pop or Multi-Pop class!", fixed = TRUE)

  pop = newPop(founderPop, simParam=SP)

  # A single Pop object is returned unchanged
  expect_identical(pop, mergeMultiPops(pop, NULL))
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

  # mergePops and mergeMultiPops do the same on MultiPop objects
  expect_identical(mergePops(mp2), mergeMultiPops(mp2))
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

  # mergePops and mergeMultiPops do the same on MultiPop objects
  expect_identical(mergePops(mp3), mergeMultiPops(mp3))
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

  # A single Pop object is returned unchanged
  expect_identical(pop, flattenMultiPop(pop))

  # Nothing to flatten, so the same input should be returned
  expect_identical(flattenMultiPop(mp1, level = 0), mp1)
  tmp = flattenMultiPop(mp1, level = 1)
  expect_identical(tmp, mp1)
  expect_identical(tmp[[1]], pop[1:2])
  expect_identical(tmp[[2]], pop[3:5])
  tmp = flattenMultiPop(mp1, level = 2)
  expect_identical(tmp, mp1)
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

  # Flattening with preserveNames
  mp_named = newMultiPop(A = pop[1:2], B = pop[3:5],
                         C = newMultiPop(D = pop[6:7], E = pop[8:10]))
  expect_identical(flattenMultiPop(mp_named, level = 1, preserveNames = "auto"),
                   newMultiPop(A = pop[1:2], B = pop[3:5], D = pop[6:7], E = pop[8:10]))
  expect_identical(flattenMultiPop(mp_named, level = 1, preserveNames = "concatenate"),
                   newMultiPop(A = pop[1:2], B = pop[3:5], C_D = pop[6:7], C_E = pop[8:10]))
  expect_identical(flattenMultiPop(mp_named, level = 1, preserveNames = "force"),
                   newMultiPop(A = pop[1:2], B = pop[3:5], C_D = pop[6:7], C_E = pop[8:10]))
  expect_identical(flattenMultiPop(mp_named, level = 1, preserveNames = "none"),
                   newMultiPop(pop[1:2], pop[3:5], pop[6:7], pop[8:10]))
  expect_identical(flattenMultiPop(mp_named, level = 2, preserveNames = "none"),
                   newMultiPop(A = pop[1:2], B = pop[3:5], 
                               C = newMultiPop(pop[6:7], pop[8:10])))
                               
  mp_named = newMultiPop(A = pop[1:2], B = pop[3:5],
                         C = newMultiPop(D = pop[6:7], D = pop[8:10]))
  expect_identical(flattenMultiPop(mp_named, level = 1, preserveNames = "auto"),
                    newMultiPop(pop[1:2], pop[3:5], pop[6:7], pop[8:10]))
  expect_identical(flattenMultiPop(mp_named, level = 1, preserveNames = "concatenate"),
                    newMultiPop(pop[1:2], pop[3:5], pop[6:7], pop[8:10]))

  expect_identical(
    expect_warning(flattenMultiPop(mp_named, level = 1, preserveNames = "force"),
                   "Duplicate names found in 'force' mode. Making names unique by appending suffixes.",
                   fixed = TRUE),
    newMultiPop(A = pop[1:2], B = pop[3:5], C_D = pop[6:7], C_D.1 = pop[8:10]))
  expect_identical(flattenMultiPop(mp_named, level = 1, preserveNames = "none"),
                    newMultiPop(pop[1:2], pop[3:5], pop[6:7], pop[8:10]))
  expect_identical(flattenMultiPop(mp_named, level = 2, preserveNames = "none"),
                    newMultiPop(A = pop[1:2], B = pop[3:5], 
                                C = newMultiPop(pop[6:7], pop[8:10])))
})

test_that("splitPop", {
  # Create founder haplotypes and set simulation parameters
  founderPop = quickHaplo(nInd = 12, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10)
  #  Create population
  pop = newPop(founderPop, simParam = SP)

  # splitPop basic behavior
  by1 = sample(LETTERS[1:3], nInd(pop), replace = TRUE)
  mp1 = splitPop(pop, by = by1)

  expect_true(isMultiPop(mp1))
  expect_length(mp1@pops, length(unique(by1)))
  expect_identical(.depthMultiPop(mp1), 1L)

  expect_identical(
    splitPop(pop, by = "A"),
    splitPop(pop, by = rep("A", length(pop)))
  )

  # splitPop recursive behavior
  mp2 = splitPop(
    pop,
    by = list(
      sample(LETTERS[1:2], nInd(pop), replace = TRUE),
      function(x) getFam(x, famType = "B")
    )
  )

  expect_true(isMultiPop(mp2))
  expect_true(all(vapply(mp2@pops, isMultiPop, logical(1))))
  expect_identical(.depthMultiPop(mp2), 2L)

  # Error handling
  expect_error(
    splitPop(newEmptyPop(ploidy = 2L, simParam = SP), by = list()),
    "`by` must have at least one grouping spec.",
    fixed = TRUE
  )

  expect_error(
    splitPop(double(), by = 1:5),
    "`x` must be a Pop or MultiPop object",
    fixed = TRUE
  )

  expect_warning(
    splitPop(pop, by = numeric(7)),
    "data length is not a multiple of split variable",
    fixed = TRUE
  )

  expect_error(
    splitPop(pop, by = list(c(NA_real_, 0))),
    "Grouping vector contains NA values.",
    fixed = TRUE
  )

  # split.default() recycling behavior
  expect_identical(
    splitPop(pop, by = 1:4),
    splitPop(pop, by = rep(1:4, length.out = nInd(pop)))
  )

  expect_error(
    splitPop(pop, by = list(list(1, 2, 3))),
    "Grouping spec must be an atomic vector or a function returning one.",
    fixed = TRUE
  )
  
  # splitPop works for MultiPop objects
  mp3 = splitPop(mp2, by = function(x) sample(LETTERS[1:2], length(x), replace = TRUE))

  expect_true(is(mp1, "MultiPop"))
  expect_true(all(vapply(mp3@pops, isMultiPop, logical(1))))
  expect_identical(.depthMultiPop(mp3), 3L)

    # Level control tests (deterministic splitter)
  mp_nested = newMultiPop(pop[1:3], newMultiPop(pop[4:6], pop[7:9]))
  splitter = function(p) rep(c("A", "B"), length.out = nInd(p))

  # level = 1: only top-level Pop nodes are split
  res1 = splitPop(mp_nested, by = splitter, level = 1)
  expect_identical(res1[[1]], splitPop(mp_nested[[1]], by = splitter))
  expect_identical(res1[[2]], mp_nested[[2]])

  # level = 2: only Pop nodes at depth 2 are split
  res2 = splitPop(mp_nested, by = splitter, level = 2)
  expect_identical(res2[[1]], mp_nested[[1]])
  expect_identical(res2[[2]], splitPop(mp_nested[[2]], by = splitter))

  # level = 1:2: Pop nodes at depth 1 and 2 are split
  res12 = splitPop(mp_nested, by = splitter, level = 1:2)
  expect_identical(res12[[1]], splitPop(mp_nested[[1]], by = splitter))
  expect_identical(res12[[2]], splitPop(mp_nested[[2]], by = splitter))

  # level = Inf: split every Pop node encountered
  resInf = splitPop(mp_nested, by = splitter, level = Inf)
  expect_identical(resInf[[1]], splitPop(mp_nested[[1]], by = splitter))
  expect_identical(resInf[[2]], splitPop(mp_nested[[2]], by = splitter))

  # Error handling for level argument
  expect_error(
    splitPop(mp_nested, by = splitter, level = -1L),
    "`level` must be a positive integer or Inf",
    fixed = TRUE
  )

  expect_error(
    splitPop(mp_nested, by = splitter, level = 4L),
    "requested `level` exceeds max depth of `x` (2)",
    fixed = TRUE
  )

  expect_error(
    splitPop(mp_nested, by = splitter, level = c(1, NA)),
    "`level` must be a numeric vector of positive integers (no NA)",
    fixed = TRUE
  )

  expect_error(
    splitPop(mp_nested, by = splitter, level = c(1, Inf)),
    "cannot mix Inf with integer levels",
    fixed = TRUE
  )

  expect_error(
    splitPop(mp_nested, by = splitter, level = c(1, 1.5)),
    "`level` must be a numeric vector of positive integers (no NA)",
    fixed = TRUE
  )

  expect_error(
    splitPop(mp_nested, by = splitter, level = c(1, 5)),
    "requested level(s) exceed max depth of `x` (2)",
    fixed = TRUE
  )

  # splitPop using a data frame for grouping
  by1 = data.frame(level1 = by1)
  expect_warning(
    expect_identical(mp1, splitPop(pop, by = by1)),
    "Mapping rows of data frame (`by`) to individuals' identifiers (`x@id`) by order.\nConsider setting row names of `by` to match `x@id` for clarity.",
    fixed = TRUE
  )
  
  rownames(by1) = letters[pop@iid]
  expect_warning(
    splitPop(pop, by = data.frame(level1 = by1)),
    "Row names of data frame `by` don't match `x@id`. Mapping rows by order instead of names.",
    fixed = TRUE
  )

  # Specify row names of the data frame to match `x@id`
  by_df = data.frame(level1 = by1, row.names = pop@id)
  expect_identical(mp1, splitPop(pop, by = by_df))
  
  # splitPop can create nested MultiPop objects using a data frame with multiple columns
  nms = sapply(mp2@pops, lengths)
  by_df = data.frame(level1 = rep(names(nms), nms),
                     level2 = "0_0",
                     row.names = mergePops(mp2)@id)
  expect_identical(mp2, splitPop(pop, by = by_df))

  # A data frame with NA values in the grouping columns creates an uneven nested structure
  by_df = data.frame(level1 = c(rep("A", 4), rep("B", 8)),
                     level2 = c(rep(NA_character_, 4), rep(LETTERS[3:4], each = 4)),
                     level3 = c(rep(NA_character_, 8), rep("E", 4)),
                     row.names = pop@id)
                     
  mp4 = newMultiPop(A = pop[1:4],
                    B = newMultiPop(
                      C = pop[5:8],
                      D = newMultiPop(E = pop[9:12])))
  
  expect_identical(mp4, splitPop(pop, by = by_df))
  

})
