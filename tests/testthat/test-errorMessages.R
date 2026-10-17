context("errorMessages")

# Input checks name the argument at fault rather than repeating the failed
# condition. These tests cover a representative set of them, one or two
# from each kind of function.

errSP = function(){
  founderPop = quickHaplo(nInd=10, nChr=2, segSites=20)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  list(founderPop=founderPop, SP=SP)
}

test_that("SimParam setup checks say what is wrong", {
  d = errSP()
  expect_error(SimParam$new(matrix(0)), "founderPop must be a MapPop")
  expect_error(d$SP$setTrackPed("yes"), "isTrackPed must be TRUE or FALSE")
  expect_error(d$SP$addTraitA(5, mean=c(0,0), var=1),
               "mean and var must have the same length")
  expect_error(d$SP$addTraitA(5, name=c("a","b")),
               "name must have one value per trait")
  expect_error(d$SP$addTraitAE(c(3,2), relAA=c(0.5,0.5)),
               "relAA must have one value per trait")
  expect_error(d$SP$addTraitAE(c(3,2)),
               "the total number of QTL must be even")
  expect_error(d$SP$setRecombRatio(0), "femaleRatio must be greater than 0")
})

test_that("missing marker names are listed", {
  d = errSP()
  expect_error(d$SP$addSnpChipByName(c("1_1","nope1","nope2")),
               "markers not found on the genetic map: nope1, nope2")
  expect_error(d$SP$importTrait(markerNames=paste0("x", 1:7),
                                addEff=rep(1,7)),
               "markerNames not found on the genetic map: x1, x2, x3, x4, x5 and 2 more")
  expect_error(d$SP$importTrait(markerNames=c("1_1","1_2"),
                                addEff=1),
               "addEff must have one row per marker in markerNames")
})

test_that("variance checks say which argument has the wrong length", {
  d = errSP()
  d$SP$addTraitA(5)
  expect_error(d$SP$setVarE(h2=c(0.5,0.5)),
               "h2 must have one value per trait")
  expect_error(d$SP$setVarE(varE=matrix(1, 2, 2)),
               "a varE matrix must have one row and one column per trait")
  expect_error(d$SP$removeTrait(3), "traits must be between 1 and nTraits")
})

test_that("import checks describe the problem with the data", {
  genMap = data.frame(id=letters[1:4], chr=c(1,1,2,2),
                      pos=c(0,0.5,0,0.5))
  geno = matrix(c(0,2,3,0, 2,0,0,2), nrow=2, byrow=TRUE,
                dimnames=list(NULL, letters[1:4]))
  expect_error(importInbredGeno(geno, genMap),
               "geno must be coded 0,1,2 or -1,0,1, but has values above 2")
  expect_error(importInbredGeno(geno*0, genMap, ped=c("a","a")),
               "ped must not contain duplicated ids")

  haplo = matrix(c(0,1,2,0, 1,0,0,1), nrow=2, byrow=TRUE,
                 dimnames=list(NULL, letters[1:4]))
  expect_error(importHaplo(haplo, genMap, ploidy=1L),
               "haplo must contain only 0 and 1")

  founderPop = quickHaplo(nInd=2, nChr=1, segSites=5)
  expect_error(sampleHaplo(founderPop, nInd=10, replace=FALSE),
               "too few haplotypes to sample without replacement")
})

test_that("population function checks name the argument", {
  d = errSP()
  d$SP$addTraitA(5)
  d$SP$addTraitA(5)
  d$SP$addSnpChip(5)
  pop = newPop(d$founderPop, simParam=d$SP)
  expect_error(editGenome(pop, ind=1, chr=1, segSites=1, allele=2,
                          simParam=d$SP),
               "allele must be 0 or 1")
  expect_error(editGenome(pop, ind=11, chr=1, segSites=1, allele=1,
                          simParam=d$SP),
               "ind must index individuals in pop")
  expect_error(setPhenoGCA(pop, testers=list(), simParam=d$SP),
               "testers must be a Pop")
  expect_error(RRBLUP2(pop, traits=1:2, use="gv", simParam=d$SP),
               "this model fits a single trait")
})
