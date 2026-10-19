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

test_that("calcPopSummary calculates summaries at different levels", {
  founderPop = quickHaplo(nInd = 32, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10, mean = c(0, 0), var = c(1, 1))

  pop = newPop(founderPop, simParam = SP)

  mp = splitPop(
    pop,
    by = data.frame(
      level1 = rep(c("G1", "G2"), each = 16),
      level2 = rep(LETTERS[1:4], each = 8),
      level3 = rep(paste0("R", 1:2), length.out = 32),
      row.names = pop@id
    )
  )

  countInd = function(x) {
    c(nInd = nInd(x))
  }

  # One summary per terminal Pop
  level3 = calcPopSummary(
    mp,
    withinPopFUN = countInd,
    acrossPopFUN = colSums,
    level = 3
  )

  expect_true(is.matrix(level3))
  expect_equal(dim(level3), c(8L, 1L))
  expect_equivalent(level3, rep(4, 8))
  expect_equal(colnames(level3), "nInd")

  source3 = attr(level3, "source")
  expect_equal(dim(source3), c(8L, 3L))
  expect_equal(
    colnames(source3),
    c("level1", "level2", "level3")
  )

  # Aggregate terminal Pops within level-2 MultiPops
  level2 = calcPopSummary(
    mp,
    withinPopFUN = countInd,
    acrossPopFUN = colSums,
    level = 2
  )

  expect_equal(dim(level2), c(4L, 1L))
  expect_equivalent(level2, rep(8, 4))
  expect_equal(dim(attr(level2, "source")), c(4L, 2L))

  # Aggregate level-2 summaries within level-1 MultiPops
  level1 = calcPopSummary(
    mp,
    withinPopFUN = countInd,
    acrossPopFUN = colSums,
    level = 1
  )

  expect_equal(dim(level1), c(2L, 1L))
  expect_equivalent(level1, c(16, 16))
  expect_equal(
    attr(level1, "source")[, 1],
    c("G1", "G2")
  )

  # Aggregate the complete hierarchy
  level0 = calcPopSummary(
    mp,
    withinPopFUN = countInd,
    acrossPopFUN = colSums,
    level = 0
  )

  expect_equal(level0, c(nInd = 32))
  expect_null(attr(level0, "source"))
})

test_that("calcPopSummary handles Pop input and function arguments", {
  founderPop = quickHaplo(nInd = 12, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)

  summary = calcPopSummary(
    pop,
    withinPopFUN = function(x, multiplier, offset) {
      c(value = multiplier * nInd(x) + offset)
    },
    withinPopArgs = list(
      multiplier = 2,
      offset = 1
    ),
    acrossPopFUN = colSums,
    level = 4
  )

  # A Pop is always summarized directly, irrespective of level
  expect_equal(summary, c(value = 25))
  expect_null(attr(summary, "source"))
})

test_that("calcPopSummary supports multiple summary values", {
  founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)

  mp = splitPop(
    pop,
    by = rep(c("A", "B"), each = 8)
  )

  result = calcPopSummary(
    mp,
    withinPopFUN = function(x) {
      c(
        nInd = nInd(x),
        meanID = mean(as.numeric(x@id))
      )
    },
    acrossPopFUN = colMeans,
    level = 1
  )

  expect_equal(dim(result), c(2L, 2L))
  expect_equal(
    colnames(result),
    c("nInd", "meanID")
  )
  expect_equal(result[, "nInd"], c(8, 8))
  expect_equal(result[, "meanID"], c(4.5, 12.5))
})

test_that("calcPopSummary returns terminal Pops below the maximum depth", {
  founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)

  mp = splitPop(
    pop,
    by = data.frame(
      rep(c("A", "B"), each = 8),
      rep(c("R1", "R2"), length.out = 16),
      row.names = pop@id
    )
  )

  result = calcPopSummary(
    mp,
    withinPopFUN = function(x) c(nInd = nInd(x)),
    acrossPopFUN = colSums,
    level = 10
  )

  expect_equal(nrow(result), 4L)
  expect_equivalent(result, rep(4, 4))

  source = attr(result, "source")
  expect_equal(ncol(source), 2L)
})

test_that("calcPopSummary handles unbalanced nesting", {
  founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)

  mp = newMultiPop(
    deep = splitPop(
      pop[1:8],
      by = rep(c("A", "B"), each = 4)
    ),
    shallow = pop[9:16]
  )

  result = calcPopSummary(
    mp,
    withinPopFUN = function(x) c(nInd = nInd(x)),
    acrossPopFUN = colSums,
    level = 2
  )

  expect_equal(nrow(result), 3L)
  expect_equivalent(result, c(4, 4, 8))

  source = attr(result, "source")
  expect_equal(
    source[, 1],
    c("deep", "deep", "shallow")
  )
  expect_equal(
    source[, 2],
    c("A", "B", NA)
  )
})

test_that("calcPopSummary validates public arguments", {
  founderPop = quickHaplo(nInd = 8, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)
  mp = splitPop(pop, by = rep(c("A", "B"), each = 4))

  expect_error(
    calcPopSummary(
      list(pop),
      withinPopFUN = nInd,
      acrossPopFUN = colSums
    ),
    "`x` must be a Pop, HybridPop or MultiPop object.",
    fixed = TRUE
  )

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = NULL,
      acrossPopFUN = colSums
    ),
    "`withinPopFUN` must be a function.",
    fixed = TRUE
  )

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = nInd,
      acrossPopFUN = NULL
    ),
    "`acrossPopFUN` must be a function.",
    fixed = TRUE
  )

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = nInd,
      withinPopArgs = NULL,
      acrossPopFUN = colSums
    ),
    "`withinPopArgs` must be a list.",
    fixed = TRUE
  )

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = nInd,
      acrossPopFUN = colSums,
      acrossPopArgs = NULL
    ),
    "`acrossPopArgs` must be a list.",
    fixed = TRUE
  )

  invalidLevels = list(
    -1,
    0.5,
    NA_real_,
    Inf,
    c(0, 1),
    "1"
  )

  for (invalidLevel in invalidLevels) {
    expect_error(
      calcPopSummary(
        mp,
        withinPopFUN = nInd,
        acrossPopFUN = colSums,
        level = invalidLevel
      ),
      "`level` must be a single non-negative integer.",
      fixed = TRUE
    )
  }
})

test_that("calcPopSummary validates within-population summaries", {
  founderPop = quickHaplo(nInd = 8, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)
  mp = splitPop(pop, by = rep(c("A", "B"), each = 4))

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = function(x) character(),
      acrossPopFUN = colSums
    ),
    paste0(
      "`withinPopFUN` must return a non-empty numeric vector or ",
      "one-row matrix without missing values."
    ),
    fixed = TRUE
  )

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = function(x) "value",
      acrossPopFUN = colSums
    ),
    paste0(
      "`withinPopFUN` must return a non-empty numeric vector or ",
      "one-row matrix without missing values."
    ),
    fixed = TRUE
  )

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = function(x) NA_real_,
      acrossPopFUN = colSums
    ),
    paste0(
      "`withinPopFUN` must return a non-empty numeric vector or ",
      "one-row matrix without missing values."
    ),
    fixed = TRUE
  )

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = function(x) matrix(1:4, nrow = 2),
      acrossPopFUN = colSums
    ),
    "`withinPopFUN` must return a vector or a one-row matrix.",
    fixed = TRUE
  )
})

test_that("calcPopSummary requires consistent summary columns", {
  founderPop = quickHaplo(nInd = 8, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)
  mp = splitPop(pop, by = rep(c("A", "B"), each = 4))

  inconsistentNames = function(x) {
    if (as.numeric(x@id[1]) <= 4) {
      c(first = 1, second = 2)
    } else {
      c(first = 1, other = 2)
    }
  }

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = inconsistentNames,
      acrossPopFUN = colSums,
      level = 1
    ),
    "Population summaries do not have consistent names.",
    fixed = TRUE
  )
})

test_that("calcPopSummary requires consistent numbers of columns", {
  founderPop = quickHaplo(nInd = 8, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)
  mp = splitPop(pop, by = rep(c("A", "B"), each = 4))

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = function(x) {
        if (as.numeric(x@id[1]) <= 4) {
          1
        } else {
          c(1, 2)
        }
      },
      acrossPopFUN = colSums,
      level = 1
    ),
    "Population summaries do not have consistent numbers of columns.",
    fixed = TRUE
  )
})

test_that("calcPopSummary validates across-population summaries", {
  founderPop = quickHaplo(nInd = 8, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)
  mp = splitPop(pop, by = rep(c("A", "B"), each = 4))

  within = function(x) {
    c(first = nInd(x), second = mean(as.numeric(x@id)))
  }

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = within,
      acrossPopFUN = function(x) "value",
      level = 0
    ),
    paste0(
      "`acrossPopFUN` must return a non-empty numeric vector or ",
      "one-row matrix without missing values."
    ),
    fixed = TRUE
  )

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = within,
      acrossPopFUN = function(x) c(first = NA_real_, second = 1),
      level = 0
    ),
    paste0(
      "`acrossPopFUN` must return a non-empty numeric vector or ",
      "one-row matrix without missing values."
    ),
    fixed = TRUE
  )

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = within,
      acrossPopFUN = function(x) matrix(1:4, nrow = 2),
      level = 0
    ),
    "`acrossPopFUN` must return a vector or a one-row matrix.",
    fixed = TRUE
  )

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = within,
      acrossPopFUN = function(x) c(first = mean(x[, "first"])),
      level = 0
    ),
    paste0(
      "`acrossPopFUN` must return one value for each column ",
      "produced by `withinPopFUN`."
    ),
    fixed = TRUE
  )

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = within,
      acrossPopFUN = function(x) {
        c(
          renamed1 = mean(x[, 1]),
          renamed2 = mean(x[, 2])
        )
      },
      level = 0
    ),
    paste0(
      "`acrossPopFUN` must preserve the names of the values ",
      "produced by `withinPopFUN`."
    ),
    fixed = TRUE
  )
})

test_that("calcPopSummary rejects empty MultiPop objects", {
  expect_error(
    calcPopSummary(
      newEmptyMultiPop(),
      withinPopFUN = nInd,
      acrossPopFUN = colSums
    ),
    "`x` contains no populations.",
    fixed = TRUE
  )
})

test_that("calcPopSummary rejects empty MultiPop objects", {
  founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)
  pop = setPheno(pop, varE = 4, simParam = SP)
  mp = splitPop(pop, by = rep(LETTERS[1:4], each = 4))
  mp[[3]] = newEmptyMultiPop()
  expect_error(
    selectMultiPop(mp, nMulti = 2, simParam = SP),
    "`x` contains an empty MultiPop object."
  )
  mp[[3]] = newEmptyPop(simParam = SP)
  expect_error(
    selectMultiPop(mp, nMulti = 2, simParam = SP),
    "`x` contains an empty Pop object."
  )
})

