context("selection")

test_that("selectInd_and_getResponse",{
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10)
  SP$setVarE(h2=0.5)
  pop = newPop(founderPop, simParam=SP)

  pop2 = selectInd(pop, 5, simParam=SP)
  expect_equal(pop2@id,
               pop[order(pop@pheno, decreasing=TRUE)[1:5]]@id)

  pop2b = selectInd(pop, 5, trait="Trait1", simParam=SP)
  expect_equal(pop2@id,
               pop2b@id)

  squaredDeviation = function(x, optima=0) (x - optima)^2
  pop3 = selectInd(pop, 5, trait=squaredDeviation, selectTop=TRUE, simParam=SP)
  expect_equal(pop3@id,
               pop[order(squaredDeviation(pop@pheno), decreasing=TRUE)[1:5]]@id)

  pop4 = selectInd(pop, 5, trait=squaredDeviation, selectTop=FALSE, simParam=SP)
  expect_equal(pop4@id,
               pop[order(squaredDeviation(pop@pheno), decreasing=FALSE)[1:5]]@id)

  pop@misc = list(smth=rnorm(10), smth2=rnorm(10))
  useFunc = function(pop, trait=NULL) pop@misc$smth + pop@misc$smth2
  pop5 = selectInd(pop, 5, use=useFunc, simParam=SP)
  expect_equal(pop5@id,
               pop[order(useFunc(pop), decreasing=TRUE)[1:5]]@id)

  useFunc2 = function(pop, trait=NULL) cbind(pop@misc$smth, pop@misc$smth2)
  trtFunc = function(x) rowSums(x)
  pop6 = selectInd(pop, 5, trait=trtFunc, use=useFunc2, simParam=SP)
  expect_equal(pop5@id, pop6@id)
})

test_that("selectPop_and_calcPopValue",{
  founderPop = quickHaplo(nInd=100, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10, mean = c(0, 0), var = c(1, 1))
  pop = newPop(founderPop, simParam=SP)
  
  # Selecting from a Pop object should return the same Pop object
  expect_equal(pop, selectPop(pop, nPop = 1, simParam = SP))

  # Multipop with 1 level of nesting
  mp1 = newMultiPop(pop[1:5], pop[6:10])
  
  # pop@pheno is empty
  expect_error(selectPop(mp1, nPop = 2, simParam = SP),
               "selection trait has missing values, phenotype may need to be set",
               fixed=TRUE)
  
  # pop@pheno is empty in the 2nd pop
  mp1@pops[[1]] = setPheno(mp1@pops[[1]], varE = c(1,1), simParam = SP)
  expect_error(selectPop(mp1, nPop = 2, simParam = SP),
               "selection trait has missing values, phenotype may need to be set",
               fixed=TRUE)
  
  
  # setPheno for all traits and pops
  mp1 = setPheno(mp1, varE = c(1,1), simParam = SP)
  
  # Suitable candidates smaller than nPop
  expect_warning(selectPop(mp1, nPop = 10, simParam = SP),
                 paste("Suitable candidate populations smaller than nPop, returning", 
                       length(mp1), "populations"))
  
  # TODO: Update support for use='bv' with genParamPop()
  # bv is not currently supported
  expect_error(selectPop(mp1, nPop = 2, use = 'bv', simParam = SP),
               "use='bv' is not currently supported for populations",
               fixed=TRUE)
  
  # Selecting nested populations from a non-nested object
  expect_error(selectPop(mp1, nPop = 1, level = 3, use = 'pheno', simParam = SP),
               paste("The MultiPop object does not contain other MultiPop objects",
                     "at this level. You may want to decrease the value of 'level'"),
               fixed=TRUE)
  
  # Use a custom trait obtained by summing up the two available traits
  values = sapply(mp1@pops, function(pop){
    mean(pop@pheno[,1:2])
  })
  expect_identical(
    selectPop(mp1, nPop = 1, simParam = SP, trait = 1:2,
              use = function(pop, trait = trait) mean(pop@pheno[,trait]))@pops[[1]],
    mp1@pops[[which.max(values)]]
  )
  expect_identical(
    selectPop(mp1, nPop = 1, selectTop = F, simParam = SP, trait = 1:2,
              use = function(pop, trait = trait) mean(pop@pheno[,trait]))@pops[[1]],
    mp1@pops[[which.min(values)]]
  )
  
  # Use a custom trait and use functions
  values = lapply(mp1@pops, function(pop, trait = 1:2){
    c(pop@pheno[,trait])
  })
  values = do.call('rbind', values)
  response = apply(values, 1, sd)
  
  expect_identical(
    selectPop(mp1, nPop = 1, simParam = SP, selectTop = F,
              use = function(pop, trait = 1:2) c(pop@pheno[,trait]),
              trait = sd)@pops[[1]],
    mp1@pops[[which.min(response)]]
  )
  expect_identical(
    selectPop(x = mp1, nPop = 1, simParam = SP, selectTop = T,
              use = function(pop, trait = 1:2) c(pop@pheno[,trait]),
              trait = sd)@pops[[1]],
    mp1@pops[[which.max(response)]]
  )
  
  # MultiPop with 1 nested object
  
  mp2 = newMultiPop(pop[1:20], pop[21:40],
                    newMultiPop(pop[41:60],
                                newMultiPop(pop[61:80], pop[81:100])))
  
  # setPheno for all traits
  mp2 = setPheno(mp2, varE = c(1,1), simParam = SP)
  
  # Selection can only be performed when all populations are Pop-class
  expect_error(
    selectPop(mp2, nPop = 1, simParam = SP),
    paste(
      "This level contains", sum(sapply(mp2@pops, isMultiPop)), 
      "MultiPop-class objects.",
      "\nSelection can only be performed when all populations at this level are",
      "Pop-class objects.\nYou may want to increase the value of 'level'"
    ), fixed = TRUE)
  expect_error(
    selectPop(mp2, nPop = 1, level = 2, simParam = SP),
    paste(
      "This level contains", sum(sapply(mp2[[3]]@pops, isMultiPop)),
      "MultiPop-class objects.",
      "\nSelection can only be performed when all populations at this level are",
      "Pop-class objects.\nYou may want to increase the value of 'level'"
    ), fixed = TRUE)
  
  # MultiPop with >1 nested object
  
  mp3 = newMultiPop(pop[1:20],
                    newMultiPop(pop[21:30],
                                newMultiPop(pop[31:35], pop[36:40])),
                    newMultiPop(pop[41:60],
                                newMultiPop(pop[61:80], pop[81:100])))
  
  mp3 = setPheno(mp3, varE = c(1,1), simParam = SP)
  
  # Selection can only be performed when all populations are Pop-class
  expect_error(
    selectPop(mp3, nPop = 1, simParam = SP),
    paste(
      "This level contains", sum(sapply(mp3@pops, isMultiPop)), 
      "MultiPop-class objects.",
      "\nSelection can only be performed when all populations at this level are",
      "Pop-class objects.\nYou may want to increase the value of 'level'"
    ), fixed = TRUE)
  expect_error(
    selectPop(mp3, nPop = 1, level = 2, simParam = SP),
    paste(
      "This level contains", sum(sapply(mp3[[3]]@pops, isMultiPop)), 
      "MultiPop-class objects.",
      "\nSelection can only be performed when all populations at this level are",
      "Pop-class objects.\nYou may want to increase the value of 'level'"
    ), fixed = TRUE)
  
  # Selection should give the same result as if we had directly selected from the appropriate nested population
  expect_identical(
    selectPop(mp3, nPop = 1, level = 3, simParam = SP)[[2]][[2]],
    selectPop(mp3[[2]][[2]], nPop = 1, level = 1, simParam = SP)
  )
  expect_identical(
    selectPop(mp3, nPop = 1, level = 3, simParam = SP)[[3]][[2]],
    selectPop(mp3[[3]][[2]], nPop = 1, level = 1, simParam = SP)
  )
})

# ---------------------------------------------------------------------------
# PART 1  Empty populations
# ---------------------------------------------------------------------------

shortfallSP = function(sexes="no"){
  set.seed(2811)
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(5)
  SP$setVarE(h2=0.5)
  SP$addSnpChip(10)
  SP$setSexes(sexes)
  pop = newPop(founderPop, simParam=SP)
  list(SP=SP, pop=pop, fam=randCross(pop, nCrosses=3, nProgeny=2,
                                     simParam=SP))
}

test_that("EMPTY a population with no individuals gives a clear error", {
  d = shortfallSP()
  empty = d$pop[integer(0)]
  msg = "pop contains no individuals to select from"
  expect_error(selectInd(empty, 1, simParam=d$SP), msg)
  expect_error(selectFam(empty, 1, simParam=d$SP), msg)
  expect_error(selectWithinFam(empty, 1, simParam=d$SP), msg)
  expect_error(selectOP(empty, 1, nSeeds=2, simParam=d$SP), msg)
  expect_error(selectOCS(empty, nInd=1, simParam=d$SP), msg)
  expect_error(selectCross(empty, nInd=1, nCrosses=1, simParam=d$SP), msg)
})

test_that("EMPTY members of a MultiPop are passed through", {
  d = shortfallSP()
  mp = newMultiPop(d$pop, d$pop[integer(0)])
  for(f in list(function(x) selectInd(x, 2, simParam=d$SP),
                function(x) selectFam(x, 1, simParam=d$SP),
                function(x) selectWithinFam(x, 1, simParam=d$SP))){
    out = f(mp)
    expect_equal(nInd(out[[2]]), 0L)
    expect_gt(nInd(out[[1]]), 0L)
  }
  out = selectInd(mp, 2, simParam=d$SP)
  expect_identical(out[[1]]@id, selectInd(d$pop, 2, simParam=d$SP)@id)
})

# ---------------------------------------------------------------------------
# PART 2  Fewer candidates than requested
# ---------------------------------------------------------------------------

test_that("SHORTFALL warning is the default and the result is unchanged", {
  d = shortfallSP()
  expect_warning(warned <- selectInd(d$pop, 20, simParam=d$SP),
                 "Suitable candidates smaller than nInd")
  quiet = selectInd(d$pop, 20, simParam=d$SP, onShortfall="none")
  expect_identical(warned@id, quiet@id)
  expect_equal(nInd(quiet), 10L)

  expect_warning(warned <- selectFam(d$fam, 5, simParam=d$SP),
                 "Suitable families smaller than nFam")
  expect_identical(warned@id,
                   selectFam(d$fam, 5, simParam=d$SP, onShortfall="none")@id)
})

test_that("SHORTFALL error stops and none is silent in every function", {
  d = shortfallSP()
  calls = list(
    function(o) selectInd(d$pop, 20, simParam=d$SP, onShortfall=o),
    function(o) selectFam(d$fam, 5, simParam=d$SP, onShortfall=o),
    function(o) selectWithinFam(d$fam, 5, simParam=d$SP, onShortfall=o),
    function(o) selectOP(d$pop, 20, nSeeds=1, simParam=d$SP,
                         onShortfall=o),
    function(o) selectOCS(d$pop, nInd=20, inbrTarget=0.5, simParam=d$SP,
                          onShortfall=o),
    function(o) selectCross(d$pop, nInd=20, nCrosses=2, simParam=d$SP,
                            onShortfall=o)
  )
  for(f in calls){
    expect_error(f("error"), "^Only [0-9]+ suitable|fewer suitable candidates than nInd = ")
    expect_warning(suppressMessages(f("none")), NA)
  }
  expect_error(selectInd(d$pop, 20, simParam=d$SP, onShortfall="error"),
               "Only 10 suitable candidates for nInd = 20")
  expect_error(selectInd(d$pop, 2, simParam=d$SP, onShortfall="stop"),
               "should be one of")
})

test_that("SHORTFALL selectWithinFam reports how many families fall short", {
  d = shortfallSP()
  expect_warning(out <- selectWithinFam(d$fam, 5, simParam=d$SP),
                 "3 of 3 families have fewer suitable candidates than nInd")
  expect_equal(nInd(out), nInd(d$fam))
  expect_silent(selectWithinFam(d$fam, 1, simParam=d$SP))
})

test_that("SHORTFALL selectFam with no eligible candidates returns none", {
  d = shortfallSP()
  expect_warning(out <- selectFam(d$fam, 1, candidates=integer(0),
                                  simParam=d$SP),
                 "returning 0 families")
  expect_equal(nInd(out), 0L)
  expect_identical(selectFam(d$fam, 1, candidates=integer(0),
                             returnPop=FALSE, simParam=d$SP,
                             onShortfall="none"),
                   integer(0))
})

test_that("SHORTFALL selectOCS stops when a sex has no candidates", {
  d = shortfallSP(sexes="yes_sys")
  females = selectInd(d$pop, 5, sex="F", simParam=d$SP)
  expect_error(selectOCS(females, nFemale=2, nMale=2, inbrTarget=0.5,
                         simParam=d$SP, onShortfall="none"),
               "no suitable male candidates")
})
