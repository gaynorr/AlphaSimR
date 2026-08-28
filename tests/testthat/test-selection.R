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

test_that("selectPop",{
  founderPop = quickHaplo(nInd=14, nChr=1, segSites=10)
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
  
  # Invalid level
  expect_error(selectPop(mp1, nPop = 2, level = 0.5, simParam = SP),
               "`level` must be a positive integer",
               fixed=TRUE)
               expect_error(selectPop(mp1, nPop = 2, level = 1:2, simParam = SP),
               "`level` must be a positive integer",
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
  
  mp2 = newMultiPop(pop[1:2], pop[3:4],
                    newMultiPop(pop[5:6],
                                newMultiPop(pop[7:8], pop[9:10])))
  
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
  
  mp3 = newMultiPop(pop[1:2],
                    newMultiPop(pop[3:4],
                                newMultiPop(pop[5:6], pop[7:8])),
                    newMultiPop(pop[9:10],
                                newMultiPop(pop[11:12], pop[13:14])))
  
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
      function(x) getFam(x, famType = "B")
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
    calcPopValue(mp1, FUN = \(x) list(x@pheno), simplify = TRUE),
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
    "`x` must be a Pop or MultiPop object.",
    fixed = TRUE
  )
  expect_error(
    calcPopValue(1:5, FUN = pheno, simplify = FALSE),
    "`x` must be a Pop or MultiPop object.",
    fixed = TRUE
  )
  expect_error(
    calcPopValue(list(pop), FUN = pheno, simplify = FALSE),
    "`x` must be a Pop or MultiPop object.",
    fixed = TRUE
  )
})
