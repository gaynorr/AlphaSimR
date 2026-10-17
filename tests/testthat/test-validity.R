context("validity")

# Validity of a Pop is checked in two places rather than one. The validity
# method holds the structural checks, which compare lengths and dimensions and
# so cost the same whatever the size of the population. R runs it on every
# object built with new(), which is every cross, every doubled haploid and
# every selection.
#
# validPopContent holds the checks that read every individual. At present that
# means the identifiers, which only fail when the names came from outside
# AlphaSimR, so it runs at newPop and not on populations derived from one.
#
# These tests exist because moving a check out of the always-run path is only
# safe if something else guarantees it still happens. They build a population
# every way the package offers and hold each one to both sets of checks.

# Returns TRUE only if the validity method finds nothing wrong
isStructurallyValid = function(object){
  isTRUE(validObject(object, test=TRUE))
}

buildSimParam = function(founderPop, sexes="no"){
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  if(sexes!="no"){
    SP$setSexes(sexes)
  }
  SP$addTraitA(10)
  SP$setVarE(h2=0.5)
  return(SP)
}

test_that("populations built by every route pass both sets of checks", {
  founderPop = quickHaplo(nInd=10, nChr=2, segSites=20)
  SP = buildSimParam(founderPop)

  pops = list()
  pops$newPop = newPop(founderPop, simParam=SP)
  base = pops$newPop

  pops$empty = newEmptyPop(simParam=SP)
  pops$setPheno = setPheno(base, simParam=SP)
  pops$randCross = randCross(base, nCrosses=5, nProgeny=2, simParam=SP)
  pops$makeCross = makeCross(base, crossPlan=cbind(1:5, 6:10), simParam=SP)
  pops$self = self(base, nProgeny=2, simParam=SP)
  pops$makeDH = makeDH(base, nDH=2, simParam=SP)
  pops$select = selectInd(base, nInd=4, simParam=SP)
  pops$subset = base[2:5]
  pops$combine = c(base[1:4], base[5:10])
  pops$merge = mergePops(list(base[1:4], base[5:10]))
  # A second generation, so that names built from earlier names are covered
  pops$secondGen = randCross(pops$randCross, nCrosses=4, simParam=SP)

  for(name in names(pops)){
    expect_true(isStructurallyValid(pops[[name]]),
                info=paste("structural checks failed for", name))
    expect_true(validPopContent(pops[[name]]),
                info=paste("content checks failed for", name))
  }

  # The whole point of the split is that an empty population is cheap and a
  # large one is not, so check that neither is treated as a special case
  expect_equal(nInd(pops$empty), 0L)
})

test_that("populations with separate sexes pass both sets of checks", {
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
  SP = buildSimParam(founderPop, sexes="yes_sys")

  pop = newPop(founderPop, simParam=SP)
  cross = randCross(pop, nCrosses=5, simParam=SP)

  expect_true(isStructurallyValid(pop))
  expect_true(validPopContent(pop))
  expect_true(isStructurallyValid(cross))
  expect_true(validPopContent(cross))
})

test_that("imported identifiers are checked and carried through", {
  haplo = rbind(c(1,1,0,1,0),
                c(1,1,0,1,0),
                c(0,1,1,0,0),
                c(0,1,1,0,0))
  colnames(haplo) = letters[1:5]
  genMap = data.frame(markerName=letters[1:5],
                      chromosome=c(1,1,1,2,2),
                      position=c(0,0.5,1,0.15,0.4))
  ped = data.frame(id=c("a","b"), mother=c(0,0), father=c(0,0))

  founderPop = importHaplo(haplo=haplo, genMap=genMap, ploidy=2L, ped=ped)
  # No traits, because this is about identifiers and two inbred individuals
  # are not enough to ask for a heritability
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L

  pop = newPop(founderPop, simParam=SP)
  expect_true(isStructurallyValid(pop))
  expect_true(validPopContent(pop))

  # Progeny take their parents' names, so the check has to hold for them too
  progeny = self(pop, nProgeny=2, simParam=SP)
  expect_true(isStructurallyValid(progeny))
  expect_true(validPopContent(progeny))

  # A name with a space is refused where the population is imported, which is
  # before any of it reaches a Pop
  badPed = data.frame(id=c("a b","c"), mother=c(0,0), father=c(0,0))
  expect_error(importHaplo(haplo=haplo, genMap=genMap, ploidy=2L,
                           ped=badPed),
               "spaces")
})

test_that("identifiers with spaces are still refused", {
  founderPop = quickHaplo(nInd=4, nChr=1, segSites=10)
  SP = buildSimParam(founderPop)
  pop = newPop(founderPop, simParam=SP)

  # Assigning a slot does not run the validity method, so these reach
  # validPopContent the way a population built outside the package would
  badId = pop
  badId@id[1] = "a b"
  expect_error(validPopContent(badId), "id can not contain spaces")

  badMother = pop
  badMother@mother[2] = "a b"
  expect_error(validPopContent(badMother), "mother can not contain spaces")

  badFather = pop
  badFather@father[3] = "a b"
  expect_error(validPopContent(badFather), "father can not contain spaces")

  # All three at once, so that the single pass over the names still reports
  # each of them
  badAll = pop
  badAll@id[1] = "a b"
  badAll@mother[1] = "c d"
  badAll@father[1] = "e f"
  msg = tryCatch(validPopContent(badAll), error=function(e) conditionMessage(e))
  expect_true(grepl("id can not contain spaces", msg, fixed=TRUE))
  expect_true(grepl("mother can not contain spaces", msg, fixed=TRUE))
  expect_true(grepl("father can not contain spaces", msg, fixed=TRUE))

})

test_that("the structural checks still catch a population that does not agree with itself", {
  founderPop = quickHaplo(nInd=6, nChr=2, segSites=10)
  SP = buildSimParam(founderPop)
  pop = newPop(founderPop, simParam=SP)

  expect_true(isStructurallyValid(pop))

  # nInd disagreeing with each of the slots it governs
  for(nm in c("sex","id","iid","mother","father","fixEff")){
    broken = pop
    slot(broken, nm, check=FALSE) = slot(broken, nm)[-1]
    expect_false(isStructurallyValid(broken),
                 info=paste("shortening", nm, "was not caught"))
  }

  # And with the matrices
  for(nm in c("gv","pheno","ebv")){
    broken = pop
    slot(broken, nm, check=FALSE) = slot(broken, nm)[-1, , drop=FALSE]
    expect_false(isStructurallyValid(broken),
                 info=paste("shortening", nm, "was not caught"))
  }

  # nTraits disagreeing with the columns of gv
  broken = pop
  broken@nTraits = 2L
  expect_false(isStructurallyValid(broken))

  # The genotypes are checked by the RawPop the Pop inherits from
  broken = pop
  broken@nChr = 1L
  expect_false(isStructurallyValid(broken))
})

test_that("MapPop and NamedMapPop validity is unchanged", {
  founderPop = quickHaplo(nInd=4, nChr=2, segSites=10)
  expect_true(isStructurallyValid(founderPop))

  haplo = rbind(c(1,1,0,1,0),
                c(1,1,0,1,0),
                c(0,1,1,0,0),
                c(0,1,1,0,0))
  colnames(haplo) = letters[1:5]
  genMap = data.frame(markerName=letters[1:5],
                      chromosome=c(1,1,1,2,2),
                      position=c(0,0.5,1,0.15,0.4))
  ped = data.frame(id=c("a","b"), mother=c(0,0), father=c(0,0))
  named = importHaplo(haplo=haplo, genMap=genMap, ploidy=2L, ped=ped)
  expect_true(isStructurallyValid(named))

  # NamedMapPop keeps its own name check, which is where imported names land
  broken = named
  broken@id[1] = "a b"
  expect_false(isStructurallyValid(broken))
})
