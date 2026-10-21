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
  
  # Invalid selectTop
  expect_error(selectPop(mp1, nPop = 2, selectTop = 1, simParam = SP),
               "`selectTop` must be TRUE or FALSE.",
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
  
  expect_error(selectPop(mp1, nPop = 2, FUN = mean, FUN.ARGS = 1, simParam = SP),
               "`FUN.ARGS` must be a list.",
               fixed=TRUE)
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
    calcPopValue(mp1, FUN = pheno, level = NA, simplify = FALSE),
    "`level` must be a single non-NA integer value.",
    fixed = TRUE
  )
  expect_error(
    calcPopValue(mp1, FUN = 1:3, simplify = FALSE),
    "`FUN` must be a function.",
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
    "`withinPopFUN` must return a numeric vector or one-row matrix.",
    fixed = TRUE
  )

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = function(x) "value",
      acrossPopFUN = colSums
    ),
    "`withinPopFUN` must return a numeric vector or one-row matrix.",
    fixed = TRUE
  )

  # Missing values are passed through rather than rejected.
  expect_equal(
    calcPopSummary(
      mp,
      withinPopFUN = function(x) c(a = NA_real_, b = 1),
      acrossPopFUN = colSums
    ),
    c(a = NA_real_, b = 2)
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

  # A 3-D array is numeric but not one row of values, for a Pop on its own
  # as well as inside a MultiPop.
  for (target in list(pop, mp)) {
    expect_error(
      calcPopSummary(
        target,
        withinPopFUN = function(x) array(1:4, c(2, 2, 1)),
        acrossPopFUN = colSums
      ),
      "`withinPopFUN` must return a vector or a one-row matrix.",
      fixed = TRUE
    )
  }
})

test_that("calcPopSummary accepts 1-D arrays as named vectors", {
  founderPop = quickHaplo(nInd = 8, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)

  # tapply() returns a 1-D array, whose names colnames() cannot read
  oneD = function(x) {
    tapply(c(nInd(x), 1), c("a", "b"), sum)
  }

  popSummary = calcPopSummary(pop, withinPopFUN = oneD, acrossPopFUN = colMeans)
  expect_null(dim(popSummary))
  expect_equal(popSummary, c(a = 8, b = 1))

  mp = newMultiPop(newMultiPop(pop[1:4], pop[5:8]), pop[1:4])
  overall = calcPopSummary(mp, withinPopFUN = oneD, acrossPopFUN = colMeans)
  expect_null(dim(overall))
  expect_equal(overall, c(a = 4, b = 1))

  byBranch = calcPopSummary(
    mp,
    withinPopFUN = oneD,
    acrossPopFUN = colMeans,
    level = 1
  )
  expect_equal(dim(byBranch), c(2L, 2L))
  expect_equal(colnames(byBranch), c("a", "b"))

  # A 1-D array from acrossPopFUN is accepted in the same way
  expect_equal(
    calcPopSummary(
      mp,
      withinPopFUN = oneD,
      acrossPopFUN = function(x) as.array(colMeans(x))
    ),
    c(a = 4, b = 1)
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
    "`acrossPopFUN` must return a numeric vector or one-row matrix.",
    fixed = TRUE
  )

  expect_equal(
    calcPopSummary(
      mp,
      withinPopFUN = within,
      acrossPopFUN = function(x) c(first = NA_real_, second = 1),
      level = 0
    ),
    c(first = NA_real_, second = 1)
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
      acrossPopFUN = function(x) array(colMeans(x), c(1, 2, 1)),
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

test_that("selectMultiPop rejects empty Pop and MultiPop members", {
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

test_that("calcPopSummary preserves scalar summary names after aggregation", {
  founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)

  mp = newMultiPop(
    newMultiPop(
      pop[1:4],
      pop[5:8],
      pop[9:12],
      pop[13:16]
    )
  )
  result = calcPopSummary(
    mp,
    withinPopFUN = function(x) {
      quantile(x@iid, probs = 0.9)
    },
    acrossPopFUN = mean,
    level = 1
  )

  expect_equal(colnames(result), "90%")

  result = calcPopSummary(
    mp,
    withinPopFUN = function(x) {
      quantile(x@iid, probs = 0.9)
    },
    acrossPopFUN = mean,
    level = 0
  )

  expect_equal(names(result), "90%")
})

test_that("calcPopSummary checks names returned by acrossPopFUN", {
  founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)

  # The nested branch is aggregated and then combined with a terminal Pop,
  # so a mislabeled aggregate would reach the overall summary.
  mp = newMultiPop(
    newMultiPop(pop[1:4], pop[5:8]),
    pop[9:16]
  )
  within = function(x) {
    c(first = nInd(x), second = mean(x@iid))
  }
  expected = calcPopSummary(
    mp,
    withinPopFUN = within,
    acrossPopFUN = colMeans,
    level = 0
  )
  expect_equal(expected, c(first = 6, second = 8.5))

  # Reordered named output is matched to the columns by name.
  reversed = calcPopSummary(
    mp,
    withinPopFUN = within,
    acrossPopFUN = function(x) rev(colMeans(x)),
    level = 0
  )
  expect_equal(reversed, expected)

  # A one-row matrix without column names takes the column names and is
  # returned as a vector.
  unnamedMatrix = calcPopSummary(
    mp,
    withinPopFUN = within,
    acrossPopFUN = function(x) matrix(colMeans(x), nrow = 1),
    level = 0
  )
  expect_null(dim(unnamedMatrix))
  expect_equal(unnamedMatrix, expected)

  # A one-row matrix from withinPopFUN is returned as a named vector.
  popSummary = calcPopSummary(
    pop,
    withinPopFUN = function(x) t(within(x)),
    acrossPopFUN = colMeans
  )
  expect_null(dim(popSummary))
  expect_equal(popSummary, within(pop))

  expect_error(
    calcPopSummary(
      mp,
      withinPopFUN = within,
      acrossPopFUN = function(x) c(a = 1, b = 2),
      level = 0
    ),
    paste0(
      "`acrossPopFUN` returned names that do not match the columns ",
      "produced by `withinPopFUN`."
    ),
    fixed = TRUE
  )

  # Empty and NA names cannot be matched by indexing with names, so
  # reordered output must still keep their values rather than give NA.
  for (secondName in c("", NA_character_)) {
    partlyNamed = function(x) {
      setNames(within(x), c("first", secondName))
    }
    expectedPartly = setNames(expected, c("first", secondName))
    expect_equal(
      calcPopSummary(
        mp,
        withinPopFUN = partlyNamed,
        acrossPopFUN = colMeans,
        level = 0
      ),
      expectedPartly
    )
    expect_equal(
      calcPopSummary(
        mp,
        withinPopFUN = partlyNamed,
        acrossPopFUN = function(x) rev(colMeans(x)),
        level = 0
      ),
      expectedPartly
    )
  }
})

test_that("selectMultiPop selects populations at different levels", {
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

  # A deterministic response based on individual IDs
  idResponse = function(pop, trait = NULL) {
    as.numeric(pop@id)
  }

  # Select the level-1 group with the highest mean ID
  selectedTop = selectMultiPop(
    mp,
    nMulti = 1,
    level = 1,
    use = idResponse,
    withinPopFUN = mean,
    acrossPopFUN = mean,
    simParam = SP
  )

  expect_equal(length(selectedTop), 1L)
  expect_equal(names(selectedTop@pops), "G2")

  # Select the level-1 group with the lowest mean ID
  selectedBottom = selectMultiPop(
    mp,
    nMulti = 1,
    level = 1,
    use = idResponse,
    selectTop = FALSE,
    withinPopFUN = mean,
    acrossPopFUN = mean,
    simParam = SP
  )

  expect_equal(length(selectedBottom), 1L)
  expect_equal(names(selectedBottom@pops), "G1")

  # Select the two highest level-2 units
  selectedLevel2 = selectMultiPop(
    mp,
    nMulti = 2,
    level = 2,
    use = idResponse,
    withinPopFUN = mean,
    acrossPopFUN = mean,
    simParam = SP
  )

  expect_equal(length(selectedLevel2), 1L)
  expect_equal(names(selectedLevel2@pops), "G2")
  expect_equal(
    names(selectedLevel2[["G2"]]@pops),
    c("C", "D")
  )
})

test_that("selectMultiPop prunes empty parent MultiPops", {
  founderPop = quickHaplo(nInd = 32, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

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

  idResponse = function(pop, trait = NULL) {
    as.numeric(pop@id)
  }

  selected = selectMultiPop(
    mp,
    nMulti = 2,
    level = 3,
    use = idResponse,
    withinPopFUN = mean,
    acrossPopFUN = mean,
    simParam = SP
  )

  # The two highest terminal responses are D/R1 and D/R2
  expect_equal(names(selected@pops), "G2")
  expect_equal(names(selected[["G2"]]@pops), "D")
  expect_equal(
    names(selected[["G2"]][["D"]]@pops),
    c("R1", "R2")
  )
})

test_that("selectMultiPop supports getResponse interfaces", {
  founderPop = quickHaplo(nInd = 32, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(
    10,
    mean = c(0, 0),
    var = c(1, 1),
    name = c("TraitA", "TraitB")
  )

  pop = newPop(founderPop, simParam = SP)
  pop = setPheno(pop, varE = c(4, 4), simParam = SP)

  mp = splitPop(
    pop,
    by = rep(c("G1", "G2", "G3", "G4"), each = 8)
  )

  byNumber = selectMultiPop(
    mp,
    nMulti = 2,
    trait = 1,
    use = "pheno",
    simParam = SP
  )

  byName = selectMultiPop(
    mp,
    nMulti = 2,
    trait = "TraitA",
    use = "pheno",
    simParam = SP
  )

  expect_equal(
    names(byNumber@pops),
    names(byName@pops)
  )

  byFunction = selectMultiPop(
    mp,
    nMulti = 2,
    trait = selIndex,
    use = "pheno",
    withinPopFUN = mean,
    acrossPopFUN = mean,
    simParam = SP,
    b = c(0.6, 0.4)
  )

  directScores = vapply(
    mp@pops,
    function(x) {
      mean(selIndex(x@pheno, b = c(0.6, 0.4)))
    },
    numeric(1)
  )

  expected = names(mp@pops)[
    order(directScores, decreasing = TRUE)[1:2]
  ]

  expect_setequal(
    names(byFunction@pops),
    expected
  )
})

test_that("selectMultiPop passes summary arguments", {
  founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)
  mp = splitPop(pop, by = data.frame(
    level1 = rep(LETTERS[1:4], each = 4),
    level2 = rep(c('R1', 'R2'), length.out = pop@nInd),
    row.names = pop@id
  ))

  idResponse = function(pop, trait = NULL) {
    as.numeric(pop@id)
  }

  selected = selectMultiPop(
    mp,
    nMulti = 1,
    use = idResponse,
    withinPopFUN = function(x, offset) mean(x) + offset,
    withinPopArgs = list(offset = 10),
    acrossPopFUN = function(x, multiplier) {
      colMeans(x) * multiplier
    },
    acrossPopArgs = list(multiplier = 2),
    simParam = SP
  )

  expect_equal(names(selected@pops), "D")
})

test_that("selectMultiPop validates selection arguments", {
  founderPop = quickHaplo(nInd = 32, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(10, mean = c(0, 0), var = c(1, 1))

  pop = newPop(founderPop, simParam = SP)
  mp = splitPop(pop, by = rep(c("G1", "G2"), each = 16))

  expect_error(
    selectMultiPop(pop, nMulti = 1, simParam = SP),
    "`x` must be a MultiPop object.",
    fixed = TRUE
  )

  invalidNMulti = list(
    0,
    -1,
    1.5,
    NA_real_,
    Inf,
    c(1, 2),
    "1"
  )

  for (value in invalidNMulti) {
    expect_error(
      selectMultiPop(mp, nMulti = value, simParam = SP),
      "`nMulti` must be a single positive integer.",
      fixed = TRUE
    )
  }

  invalidLevels = list(
    0,
    -1,
    1.5,
    NA_real_,
    Inf,
    c(1, 2),
    "1"
  )

  for (value in invalidLevels) {
    expect_error(
      selectMultiPop(
        mp,
        nMulti = 1,
        level = value,
        simParam = SP
      ),
      "`level` must be a single positive integer.",
      fixed = TRUE
    )
  }

  invalidSelectTop = list(
    1,
    NA,
    c(TRUE, FALSE),
    "TRUE"
  )

  for (value in invalidSelectTop) {
    expect_error(
      selectMultiPop(
        mp,
        nMulti = 1,
        selectTop = value,
        simParam = SP
      ),
      "`selectTop` must be TRUE or FALSE.",
      fixed = TRUE
    )
  }

  expect_error(
    selectMultiPop(
      mp,
      nMulti = 1,
      withinPopFUN = pi,
      simParam = SP
    ),
    "`withinPopFUN` must be a function.",
    fixed = TRUE
  )

  expect_error(
    selectMultiPop(
      mp,
      nMulti = 1,
      withinPopArgs = NULL,
      simParam = SP
    ),
    "`withinPopArgs` must be a list.",
    fixed = TRUE
  )

  expect_error(
    selectMultiPop(
      mp,
      nMulti = 1,
      use = "bv",
      simParam = SP
    ),
    "use='bv' is not currently supported for populations.",
    fixed = TRUE
  )

  expect_error(
    selectMultiPop(
      mp,
      nMulti = 1,
      acrossPopFUN = NULL,
      simParam = SP
    ),
    "`acrossPopFUN` must be a function.",
    fixed = TRUE
  )

  expect_error(
    selectMultiPop(
      mp,
      nMulti = 1,
      acrossPopArgs = NULL,
      simParam = SP
    ),
    "`acrossPopArgs` must be a list.",
    fixed = TRUE
  )
})

test_that("selectMultiPop requires one finite score per candidate", {
  founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)

  idResponse = function(pop, trait = NULL) {
    as.numeric(pop@id)
  }

  mp = splitPop(pop, by = rep(LETTERS[1:4], each = 4))

  expect_error(
    selectMultiPop(
      mp,
      nMulti = 1,
      use = idResponse,
      withinPopFUN = function(x) {
        c(mean = mean(x), variance = var(x))
      },
      acrossPopFUN = colMeans,
      simParam = SP
    ),
    paste0(
      "`withinPopFUN` and `acrossPopFUN` must produce exactly one ",
      "selection value per candidate."
    ),
    fixed = TRUE
  )

  expect_error(
    selectMultiPop(
      mp,
      nMulti = 1,
      use = idResponse,
      withinPopFUN = function(x) Inf,
      acrossPopFUN = mean,
      simParam = SP
    ),
    "Candidate selection values must be finite and non-missing.",
    fixed = TRUE
  )
})

test_that("selectMultiPop handles candidate counts", {
  founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)
  mp = splitPop(pop, by = rep(LETTERS[1:4], each = 4))

  idResponse = function(pop, trait = NULL) {
    as.numeric(pop@id)
  }

  expect_equal(
    expect_warning(
      selectMultiPop(mp, nMulti = 5, use = idResponse, simParam = SP),
      paste0(
        "`nMulti` exceeds the number of available candidates; ",
        "returning all candidates."
      ),
      fixed = TRUE
    ),
    mp
  )

  expect_equal(
    selectMultiPop(
      mp,
      nMulti = 4,
      use = idResponse,
      simParam = SP
    ),
    mp
  )
})

test_that("selectMultiPop resolves ties in population order", {
  founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)
  mp = splitPop(pop, by = rep(LETTERS[1:4], each = 4))

  constantResponse = function(pop, trait = NULL) {
    rep(1, nInd(pop))
  }

  selected = selectMultiPop(
    mp,
    nMulti = 2,
    use = constantResponse,
    withinPopFUN = mean,
    acrossPopFUN = mean,
    simParam = SP
  )

  expect_equal(
    names(selected@pops),
    c("A", "B")
  )
})

test_that("selectMultiPop propagates getResponse errors", {
  founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)
  mp = splitPop(pop, by = rep(LETTERS[1:4], each = 4))

  expect_error(
    selectMultiPop(
      mp,
      nMulti = 2,
      use = "notAResponse",
      simParam = SP
    ),
    "Use=notaresponse is not an option",
    fixed = TRUE
  )

  expect_error(
    selectMultiPop(
      mp,
      nMulti = 2,
      trait = "unknownTrait",
      simParam = SP
    ),
    "'unknownTrait' did not match any trait names",
    fixed = TRUE
  )

  # No phenotypes have been set
  expect_error(
    selectMultiPop(
      mp,
      nMulti = 2,
      use = "pheno",
      simParam = SP
    ),
    "selection trait has missing values, phenotype may need to be set",
    fixed = TRUE
  )
})

test_that("selectMultiPop uses positional paths for partially named objects", {
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

  # Remove one name to create a partially named hierarchy
  names(mp)[2] = NA

  useID = function(pop, trait = NULL) {
    as.numeric(pop@id)
  }

  selected = selectMultiPop(
    mp,
    nMulti = 1,
    level = 2,
    use = useID,
    selectTop = FALSE,
    simParam = SP
  )

  expect_true(isMultiPop(selected))
  expect_equal(
    sum(vapply(selected@pops, length, integer(1))),
    1L
  )
})

test_that("selectMultiPop removes multiple sibling candidates", {
  founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)

  mp = newMultiPop(
    newMultiPop(
      pop[1:4],
      pop[5:8],
      pop[9:12],
      pop[13:16]
    )
  )

  idResponse = function(pop, trait = NULL) {
    as.numeric(pop@id)
  }

  selected = selectMultiPop(
    mp,
    nMulti = 2,
    level = 2,
    use = idResponse,
    selectTop = TRUE,
    simParam = SP
  )

  expect_equal(length(selected[[1]]), 2L)
  expect_equal(selected[[1]][[1]]@id, pop[9:12]@id)
  expect_equal(selected[[1]][[2]]@id, pop[13:16]@id)
})

test_that("selectMultiPop prunes a flat hierarchy in one traversal", {
  founderPop = quickHaplo(
    nInd = 40,
    nChr = 1,
    segSites = 10
  )
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)

  mp = newMultiPop(
    A = pop[1:4],
    B = pop[5:8],
    C = pop[9:12],
    D = pop[13:16],
    E = pop[17:20],
    F = pop[21:24],
    G = pop[25:28],
    H = pop[29:32],
    I = pop[33:36],
    J = pop[37:40]
  )

  useID = function(pop, trait = NULL) {
    as.numeric(pop@id)
  }

  selected = selectMultiPop(
    mp,
    nMulti = 1,
    level = 1,
    use = useID,
    simParam = SP
  )

  expect_equal(length(selected), 1L)
  expect_equal(names(selected), "J")
  expect_equal(selected[[1]]@id, pop[37:40]@id)
})

test_that("selectMultiPop prunes nested paths in one traversal", {
  founderPop = quickHaplo(
    nInd = 32,
    nChr = 1,
    segSites = 10
  )
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)

  mp = splitPop(
    pop,
    by = data.frame(
      level1 = rep(c("G1", "G2"), each = 16),
      level2 = rep(LETTERS[1:4], each = 8),
      level3 = rep(c("R1", "R2"), length.out = 32),
      row.names = pop@id
    )
  )

  useID = function(pop, trait = NULL) {
    as.numeric(pop@id)
  }

  selected = selectMultiPop(
    mp,
    nMulti = 2,
    level = 3,
    use = useID,
    simParam = SP
  )

  expect_equal(names(selected), "G2")
  expect_equal(names(selected[["G2"]]), "D")
  expect_equal(
    names(selected[["G2"]][["D"]]),
    c("R1", "R2")
  )
})

test_that(".stripPops removes multiple paths and prunes empty parents", {
  founderPop = quickHaplo(
    nInd = 24,
    nChr = 1,
    segSites = 10
  )
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)

  mp = newMultiPop(
    first = newMultiPop(
      A = pop[1:4],
      B = pop[5:8]
    ),
    second = newMultiPop(
      C = pop[9:12],
      D = pop[13:16]
    ),
    third = newMultiPop(
      E = pop[17:20],
      F = pop[21:24]
    )
  )

  result = AlphaSimR:::.stripPops(
    mp,
    paths = rbind(
      c(1L, 1L),
      c(1L, 2L),
      c(2L, 1L)
    )
  )

  # The first parent becomes empty and is removed.
  expect_equal(
    names(result),
    c("second", "third")
  )

  # The remaining children retain their original names and order.
  expect_equal(
    names(result[["second"]]),
    "D"
  )
  expect_equal(
    names(result[["third"]]),
    c("E", "F")
  )

  expect_equal(
    result[["second"]][["D"]]@id,
    pop[13:16]@id
  )
})

test_that(".stripPops handles paths ending at different levels", {
  founderPop = quickHaplo(
    nInd = 16,
    nChr = 1,
    segSites = 10
  )
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)

  pop = newPop(founderPop, simParam = SP)

  mp = newMultiPop(
    first = newMultiPop(
      A = pop[1:4],
      B = pop[5:8]
    ),
    second = pop[9:12],
    third = pop[13:16]
  )

  result = AlphaSimR:::.stripPops(
    mp,
    paths = rbind(
      c(1L, 2L),
      c(2L, NA_integer_)
    )
  )

  expect_equal(
    names(result),
    c("first", "third")
  )
  expect_equal(
    names(result[["first"]]),
    "A"
  )
})
