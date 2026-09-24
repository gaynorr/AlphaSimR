context("misc")

test_that("misc_and_miscPop", {
  founderPop = quickHaplo(nInd = 2, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10)
  popOrig = newPop(founderPop, simParam = SP)
  multiPop = newMultiPop(popOrig, popOrig)

  expect_equal(popOrig@misc, list())
  expect_equal(popOrig@miscPop, list())

  popSub = popOrig[1]
  expect_equal(popSub@misc, list())
  expect_equal(popSub@miscPop, list())

  pop = popOrig
  pop@misc$vec = rnorm(n = 2)
  pop@misc$mat = matrix(1:4, nrow = 2)
  pop@misc$mtP = popOrig # hmm, should this actually be multiple pop objects or one with multiple individuals?
  pop@misc$mtLP = list(popOrig, popOrig)
  pop@misc$mtMP = multiPop
  # setting these miscPop elements just as an example - we are not testing them below,
  # because they get dropped in most/all operations on a pop
  pop@miscPop$vec = sum(pop@misc$vec)
  pop@miscPop$af = colMeans(pullSegSiteGeno(pop, simParam = SP))
  pop@miscPop$mat = matrix(1:4, nrow = 2)

  popSub = pop[1]
  expect_equal(popSub@misc$vec, pop@misc$vec[1])
  expect_equal(popSub@misc$mat, pop@misc$mat[1, , drop = FALSE])
  expect_equal(popSub@misc$mtP, pop@misc$mtP[1])
  expect_equal(popSub@misc$mtLP, pop@misc$mtLP[1])
  expect_equal(popSub@misc$mtMP, pop@misc$mtMP[1])
  expect_equal(popSub@miscPop, list())

  popSub = pop[0]
  expect_equal(popSub@misc$vec, numeric(0))
  expect_equal(popSub@misc$mat, pop@misc$mat[0, , drop = FALSE])
  expect_equal(popSub@misc$mtP, newEmptyPop(simParam = SP))
  expect_equal(popSub@misc$mtLP, list())
  expect_equal(popSub@misc$mtMP, new("MultiPop", pops = list()))
  expect_equal(popSub@miscPop, list())

  popC = c(pop, pop)
  expect_equal(popC@misc$vec, c(pop@misc$vec, pop@misc$vec))
  expect_equal(popC@misc$mat, rbind(pop@misc$mat, pop@misc$mat))
  expect_equal(popC@misc$mtP, c(pop@misc$mtP, pop@misc$mtP))
  expect_equal(popC@misc$mtLP, c(pop@misc$mtLP, pop@misc$mtLP))
  expect_equal(popC@misc$mtMP, c(pop@misc$mtMP, pop@misc$mtMP))
  expect_equal(popC@miscPop, list())

  popC = c(pop, pop[1])
  expect_equal(popC@misc$vec, c(pop@misc$vec, pop@misc$vec[1]))
  expect_equal(popC@misc$mat, rbind(pop@misc$mat, pop@misc$mat[1, ]))
  expect_equal(popC@misc$mtP, c(pop@misc$mtP, pop@misc$mtP[1]))
  expect_equal(popC@misc$mtLP, c(pop@misc$mtLP, pop@misc$mtLP[1]))
  expect_equal(popC@misc$mtMP, c(pop@misc$mtMP, pop@misc$mtMP[1]))
  expect_equal(popC@miscPop, list())

  popA = popB = pop
  popA@misc = pop@misc[1]
  popB@misc = pop@misc[2]
  expect_warning(
    c(popA, popB)@misc,
    regexp = "misc element names do not match - setting misc to an empty list!"
  )

  popA@misc = pop@misc[1:2]
  popB@misc = pop@misc[2:1]
  expect_warning(
    c(popA, popB)@misc,
    regexp = "misc element names do not match - setting misc to an empty list!"
  )

  popA@misc = list()
  popB@misc = pop@misc[2:1]

  popList = list(popA, popB)
  expect_warning(
    c(popA, popB)@misc,
    regexp = "number of misc elements differs - setting misc to an empty list!"
  )

  popA@misc = pop@misc[2:1]
  popB@misc = list()
  expect_warning(
    c(popA, popB)@misc,
    regexp = "number of misc elements differs - setting misc to an empty list!"
  )
})

test_that("mutate", {
  founderPop = newMapPop(
    list(c(0, 0, 0)),
    list(matrix(c(0, 0, 0, 0, 0, 0), nrow = 2, ncol = 3))
  )
  SP = SimParam$new(founderPop = founderPop)
  SP$nThreads = 1L
  pop = newPop(founderPop, simParam = SP)
  hapBefore = pullSegSiteHaplo(pop, simParam = SP)
  pop = mutate(pop, mutRate = 0, simParam = SP)
  hapAfter = pullSegSiteHaplo(pop, simParam = SP)
  expect_true(sum(hapAfter - hapBefore) == 0)
  pop = mutate(pop, mutRate = 1, simParam = SP)
  hapAfter = pullSegSiteHaplo(pop, simParam = SP)
  expect_true(sum(hapAfter - hapBefore) == 6)
})

test_that("mutate is a generic that forwards what it cannot handle", {
  expect_true(isGeneric("mutate"))
  expect_true(existsMethod("mutate", "Pop"))
  expect_true(existsMethod("mutate", "ANY"))

  # A data.frame is not a population, so it goes to the next mutate on the
  # search path rather than to an error. The call is qualified because an
  # attached mutate would otherwise be found before the generic, which is
  # exactly the situation being modelled.
  other = new.env()
  assign("mutate", function(pop, ...) "forwarded", envir = other)
  attach(other, name = "testOtherMutate", warn.conflicts = FALSE)
  on.exit(detach("testOtherMutate"), add = TRUE)

  expect_equal(AlphaSimR::mutate(data.frame(x = 1)), "forwarded")

  # A population still reaches the Pop method while the other one is
  # attached, because dispatch happens inside the generic
  founderPop = quickHaplo(nInd = 2, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  pop = newPop(founderPop, simParam = SP)
  expect_s4_class(AlphaSimR::mutate(pop, mutRate = 0, simParam = SP), "Pop")
})

test_that("mutate says so when there is nothing to forward to", {
  # Only meaningful when nothing else on the search path provides one
  hasOther = vapply(search(), function(w){
    !identical(w, "package:AlphaSimR") &&
      exists("mutate", envir = as.environment(w), inherits = FALSE)
  }, logical(1))
  skip_if(any(hasOther), "another mutate is attached")

  expect_error(AlphaSimR::mutate(data.frame(x = 1)), "expects a Pop")
})

test_that("NamedMapPop subsetting by id", {
  mapPop = quickHaplo(nInd=4, nChr=2, segSites=10)
  named = new("NamedMapPop",
              id=c("a","b","c","d"),
              mother=c("0","0","a","b"),
              father=c("0","0","b","a"),
              mapPop)

  # A single id
  expect_equal(named["c"]@id, "c")
  expect_equal(named["c"]@nInd, 1L)

  # Several ids, in the order given rather than the order stored
  sub = named[c("d","a")]
  expect_equal(sub@id, c("d","a"))
  expect_equal(sub@mother, c("b","0"))
  expect_equal(sub@father, c("a","0"))
  expect_equal(sub@nInd, 2L)

  # The genotypes follow the ids
  expect_equal(pullSegSiteGeno(sub),
               pullSegSiteGeno(named)[c("d","a"), , drop=FALSE])

  # Integer indexing is unchanged, and agrees with the ids
  expect_equal(pullSegSiteGeno(named[c(4,1)]), pullSegSiteGeno(sub))

  # A repeated id is a repeated individual, as it is for a Pop
  expect_equal(named[c("a","a")]@id, c("a","a"))

  # Unknown ids and out of range indices are both refused
  expect_error(named["z"], "invalid individuals")
  expect_error(named[c("a","z")], "invalid individuals")
  expect_error(named[5], "invalid individuals")

  # A MapPop carries no ids, so it stays index only
  expect_error(mapPop["a"])
})
