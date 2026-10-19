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

  expect_error(selectPop(mp1, nPop = 2, level = 0.5, simParam = SP),
               "`level` must be a single non-NA integer value.",
               fixed=TRUE)
  expect_error(selectPop(mp1, nPop = 2, level = 1:2, simParam = SP),
               "`level` must be a single non-NA integer value.",
               fixed=TRUE)
  
  # TODO: Update support for use='bv' with genParamPop()
  # bv is not currently supported
  expect_error(selectPop(mp1, nPop = 2, use = 'bv', simParam = SP),
               "use='bv' is not currently supported for populations",
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

test_that("SELECTION forwards functions for families and MultiPop", {
  f = miscFunctionFixture()
  pop = f$pop
  score = function(pop, ...) pop@misc$value
  fam = selectFam(pop, nFam=1, use=score, simParam=f$sp)
  expect_identical(fam@miscPop$total, sum(fam@misc$value))
  within = selectWithinFam(pop, nInd=1, use=score, simParam=f$sp)
  expect_identical(within@miscPop$total, sum(within@misc$value))
  mp = selectInd(newMultiPop(pop, pop), nInd=1, use=score, simParam=f$sp)
  expect_identical(mp[[1]]@miscPop$total, 8)
  expect_identical(mp[[2]]@miscPop$total, 8)
})

# tests/testthat/test-calcPopValue.R
test_that("calcPopValue", {
  founderPop = quickHaplo(nInd = 14, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10, mean = c(0, 0), var = c(1, 1))
  SP$setVarE(h2 = c(0.6, 0.4))
  pop = newPop(founderPop, simParam = SP)

  # Pop case
  expect_identical(
    calcPopValue(pop, FUN = pheno, simplify = FALSE),
    pheno(pop)
  )
  expect_identical(
    calcPopValue(
      pop,
      FUN = function(x, ...) cor(x@pheno, ...),
      simplify = TRUE,
      method = "kendall"
    ),
    cor(pheno(pop), method = "kendall")
  )

  # Nested MultiPop case
  mp1 = splitPop(
    randCross(pop, nCrosses = 3, nProgeny = 4, simParam = SP),
    by = list(
      function(x) rep(LETTERS[1:2], length.out = length(x)),
      function(x) AlphaSimR:::getFam(x, famType = "B")
    )
  )

  # simplify=FALSE preserves nesting and names
  outF = calcPopValue(mp1, FUN = pheno, simplify = FALSE)
  expect_true(is.list(outF))
  expect_identical(names(outF), names(mp1))
  expect_true(is.list(outF$A))
  expect_identical(outF$A, lapply(mp1$A@pops, pheno))
  expect_identical(outF$B, lapply(mp1$B@pops, pheno))

  # level > depth returns unchanged structure
  expect_identical(
    outF,
    calcPopValue(mp1, FUN = pheno, simplify = TRUE, level = 3)
  )

  # simplify=TRUE level=2 attaches source with level2
  out2 = calcPopValue(
    mp1,
    FUN = pheno,
    simplify = TRUE,
    level = 2
  )
  expect_true(is.list(out2))
  expect_identical(names(out2), names(mp1))

  outA = calcPopValue(mp1$A, FUN = pheno, simplify = TRUE)
  expect_identical(attributes(out2$A)$source[[1]], attributes(outA)$source[[1]])
  expect_equal(nrow(attributes(outA)$source), nrow(outA))

  outB = calcPopValue(mp1$B, FUN = pheno, simplify = TRUE)
  expect_identical(attributes(out2$B)$source[[1]], attributes(outB)$source[[1]])
  expect_equal(nrow(attributes(outB)$source), nrow(outB))

  expect_equal(rbind(outA), do.call('rbind', lapply(mp1$A@pops, pheno)))
  expect_equal(rbind(outB), do.call('rbind', lapply(mp1$B@pops, pheno)))

  # simplify=TRUE level=1 attaches level1/level2 source
  out1 = calcPopValue(
    mp1,
    FUN = pheno,
    simplify = TRUE,
    level = 1
  )
  expect_true(is.matrix(out1))
  expect_equal(nrow(out1), nInd(mergeMultiPops(mp1)))
  expect_identical(
    rbind(attributes(out2$A)$source, attributes(out2$B)$source),
    attributes(out1)$source[, 2, drop = FALSE]
  )
  expect_equal(rbind(out1), rbind(outA, outB))

  # Ragged and nested MultiPop
  mp3 = newMultiPop(
    p1 = pop[1:2],
    mp1 = newMultiPop(
      p2 = pop[3:4],
      mp2 = newMultiPop(p3 = pop[5:6], p4 = pop[7:8])
    ),
    mp3 = newMultiPop(
      p5 = pop[9:10],
      mp4 = newMultiPop(p6 = pop[11:12], p7 = pop[13:14])
    )
  )

  # simplify=FALSE preserves nesting and names
  outF = calcPopValue(mp3, FUN = pheno, simplify = FALSE)
  expect_true(is.list(outF))
  expect_identical(names(outF), names(mp3))
  expect_identical(names(outF$mp1), names(mp3$mp1))
  expect_identical(names(outF$mp3), names(mp3$mp3))
  expect_identical(outF$mp1$mp2, lapply(mp3$mp1$mp2@pops, pheno))
  expect_identical(outF$mp3$mp4, lapply(mp3$mp3$mp4@pops, pheno))

  # simplify=TRUE level=1 attaches level1/level2/level3 source
  out1 = calcPopValue(
    mp3,
    FUN = pheno,
    simplify = TRUE,
    level = 1
  )
  expect_true(is.matrix(out1))
  expect_equal(nrow(out1), nInd(mergeMultiPops(mp3)))

  # Unnamed multipop
  un_mp3 = newMultiPop(
    pop[1:2],
    newMultiPop(
      pop[3:4],
      newMultiPop(pop[5:6], pop[7:8])
    ),
    newMultiPop(
      pop[9:10],
      newMultiPop(pop[11:12], pop[13:14])
    )
  )

  # simplify=FALSE preserves nesting
  outF = calcPopValue(un_mp3, FUN = pheno, simplify = FALSE)
  expect_true(is.list(outF))
  expect_identical(names(outF), names(un_mp3))
  expect_identical(outF[[1]], un_mp3[[1]]@pheno)
  expect_identical(outF[[2]][[1]], un_mp3[[2]][[1]]@pheno)
  expect_identical(outF[[2]][[2]][[1]], un_mp3[[2]][[2]][[1]]@pheno)
  expect_identical(outF[[2]][[2]][[2]], un_mp3[[2]][[2]][[2]]@pheno)
  expect_identical(outF[[3]][[1]], un_mp3[[3]][[1]]@pheno)
  expect_identical(outF[[3]][[2]][[1]], un_mp3[[3]][[2]][[1]]@pheno)
  expect_identical(outF[[3]][[2]][[2]], un_mp3[[3]][[2]][[2]]@pheno)

  # simplify=TRUE level=1 attaches level1/level2/level3 source
  out1 = calcPopValue(
    un_mp3,
    FUN = pheno,
    simplify = TRUE,
    level = 1
  )
  expect_true(is.matrix(out1))
  expect_equal(nrow(out1), nInd(mergeMultiPops(un_mp3)))

  # Error and warning handling
  expect_warning(
    calcPopValue(mp1, FUN = function(x) list(x@pheno), simplify = TRUE),
    "Some values returned by FUN have unsupported types for simplification. Returning list output.",
    fixed = TRUE
  )
  colnames(mp1$A[[1]]@pheno) = c("tRAIT2","tRAIT1")
  expect_warning(
    calcPopValue(mp1, FUN = pheno, simplify = TRUE),
    "Some values returned by FUN do not have consistent column names. Returning list output.",
    fixed = TRUE
  )
  expect_error(
    calcPopValue("not_a_pop", FUN = pheno, simplify = FALSE),
    "`x` must be a Pop, HybridPop or MultiPop object.",
    fixed = TRUE
  )
  expect_error(
    calcPopValue(1:5, FUN = pheno, simplify = FALSE),
    "`x` must be a Pop, HybridPop or MultiPop object.",
    fixed = TRUE
  )
  expect_error(
    calcPopValue(list(pop), FUN = pheno, simplify = FALSE),
    "`x` must be a Pop, HybridPop or MultiPop object.",
    fixed = TRUE
  )
})

test_that("calcPopValue KEEPS names when simplification falls back to a list", {
  founderPop = quickHaplo(nInd=5, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10)
  pop = newPop(founderPop, simParam=SP)
  mp = newMultiPop(A=pop[1:2], B=pop[3:5])

  # A list is not a type that can be simplified.
  out = expect_warning(
    calcPopValue(mp, FUN=function(x) list(x@nInd), simplify=TRUE),
    "unsupported types"
  )
  expect_identical(names(out), c("A", "B"))
  expect_identical(out$A, list(2L))

  # Vectors whose names differ cannot be bound by column.
  out = expect_warning(
    calcPopValue(
      mp,
      FUN=function(x) setNames(x@nInd, if(x@nInd == 2L) "a" else "b"),
      simplify=TRUE
    ),
    "consistent column names"
  )
  expect_identical(names(out), c("A", "B"))
  expect_identical(out$B, c(b=3L))
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

test_that("UNKNOWN parents of 0 or NA each form their own family", {
  d = shortfallSP()

  # Known parents keep the labels they have always had
  expect_identical(AlphaSimR:::getFam(d$fam, "B"),
                   paste(d$fam@mother, d$fam@father, sep="_"))
  expect_identical(AlphaSimR:::getFam(d$fam, "F"), d$fam@mother)
  expect_identical(AlphaSimR:::getFam(d$fam, "M"), d$fam@father)

  # The case from discussion #71: with the parents blanked, every
  # individual used to fall into one family and only one was returned
  for(unknown in c("0", NA_character_)){
    fam = d$fam
    fam@mother[] = unknown
    fam@father[] = unknown
    for(famType in c("B", "F", "M")){
      expect_identical(AlphaSimR:::getFam(fam, famType),
                       paste0("unknown_", fam@id))
      expect_equal(nInd(selectWithinFam(fam, 1, famType=famType,
                                        simParam=d$SP)),
                   nInd(fam))
    }
  }

  # Founders have parents of "0", so each founder is a family
  expect_equal(nInd(selectFam(d$pop, 2, simParam=d$SP)), 2L)

  # pedigreeCross records unknown parents as NA. Individual 4 has a known
  # mother and an unknown father, so it is a half sib of 5 on the female
  # side but cannot be placed in a full-sib family.
  id = as.character(1:5)
  mother = c(NA, NA, NA, "1", "1")
  father = c(NA, NA, "2", NA, "2")
  ped = pedigreeCross(d$pop, id, mother, father, simParam=d$SP)
  expect_identical(AlphaSimR:::getFam(ped, "B"),
                   c(paste0("unknown_", 1:4), "1_2"))
  expect_identical(AlphaSimR:::getFam(ped, "F"),
                   c(paste0("unknown_", 1:3), "1", "1"))
  expect_identical(AlphaSimR:::getFam(ped, "M"),
                   c("unknown_1", "unknown_2", "2", "unknown_4", "2"))

  # Nothing is dropped: aggregate would discard NA labels in selectFam
  expect_equal(nInd(selectFam(ped, 5, simParam=d$SP)), 5L)
  expect_equal(nInd(selectWithinFam(ped, 1, famType="F", simParam=d$SP)),
               4L)
})


test_that("HYBRID selection preserves population indexing", {
  founder = quickHaplo(nInd=4, nChr=1, segSites=10)
  sp = SimParam$new(founder)
  sp$nThreads = 1L
  sp$addTraitA(5)
  parents = newPop(founder, simParam=sp)
  hybrid = hybridCross(parents, parents, returnHybridPop=TRUE, simParam=sp)
  # Fixed scores make selection independent of the sampled genetic values.
  score = function(pop, ...) seq_len(pop@nInd)
  for(select in list(selectInd, selectFam, selectWithinFam)){
    index = select(hybrid, 1, use=score, returnPop=FALSE, simParam=sp)
    out = select(hybrid, 1, use=score, simParam=sp)
    expect_s4_class(out, "HybridPop")
    expect_identical(out, hybrid[index])
  }
  empty = selectFam(hybrid, 1, candidates=integer(0),
                    simParam=sp, onShortfall="none")
  expect_identical(empty, hybrid[integer(0)])
})

test_that("EMPTY family selection uses explicit miscellaneous hooks", {
  f = miscFunctionFixture()
  check = function(){
    out = selectFam(f$pop, 1, candidates=integer(0),
                    simParam=f$sp, onShortfall="none")
    expect_equal(out@nInd, 0L)
    expect_identical(out@miscPop, list(total=0, ids=character(0)))
    expect_identical(out@misc$value, numeric(0))
  }
  withWrongGlobalMiscSP(f$sp, check())
  withoutGlobalSimParam(check())
})
