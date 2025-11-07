context("selection")

test_that("selectInd_and_getResponse",{
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
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

test_that("selectPop_and_getResponsePop",{
  founderPop = quickHaplo(nInd=100, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10, mean = c(0, 0), var = c(1, 1))
  SP$setVarE(h2=rep(0.5, 2))
  pop = newPop(founderPop, simParam=SP)
  
  # Multipop with 1 level of nesting
  
  multiPop = newMultiPop(pop[1:5], pop[6:10])
  
  # miscPop@pheno is NULL
  expect_error(selectPop(multiPop, nPop = 2, simParam = SP),
               "One or more populations does not have a valid pop@miscPop$pheno slot",
               fixed=TRUE)
  
  # miscPop@pheno is NULL in the 2nd pop
  multiPop = newMultiPop(setPhenoPop(pop[1:5], simParam = SP), pop[6:10])
  expect_error(selectPop(multiPop, nPop = 2, simParam = SP),
               "One or more populations does not have a valid pop@miscPop$pheno slot",
               fixed=TRUE)
  
  
  # miscPop@pheno is empty in one pop
  multiPop@pops[[2]]@miscPop$pheno = matrix(NA, ncol = 2)
  expect_error(selectPop(multiPop, nPop = 2, simParam = SP),
               "One or more populations has an emtpy pop@miscPop$pheno matrix",
               fixed=TRUE)
  
  # miscPop@pheno is empty in all pops
  multiPop@pops[[1]]@miscPop$pheno = matrix(NA, ncol = 2)
  expect_error(selectPop(multiPop, nPop = 2, simParam = SP),
               "One or more populations has an emtpy pop@miscPop$pheno matrix",
               fixed=TRUE)
  
  # setPhenoPop for all traits
  multiPop = setPhenoPop(multiPop, simParam = SP)
  
  # Suitable candidates smaller than nPop
  expect_warning(selectPop(multiPop, nPop = 10, simParam = SP),
                 paste("Suitable candidate populations smaller than nPop, returning", 
                       length(multiPop), "populations"))
  
  # gv is not present in the population
  expect_error(selectPop(multiPop, nPop = 2, use = 'gv', simParam = SP),
               "One or more populations does not have a valid pop@miscPop$gv slot",
               fixed=TRUE)
  
  # Selecting nested populations from a non-nested object
  expect_error(selectPop(multiPop, nPop = 1, level = 3, use = 'pheno', simParam = SP),
               paste("The MultiPop object does not contain other MultiPop objects",
                     "at this level. You may want to decrease the value of 'level'"),
               fixed=TRUE)
  
  # Use a custom trait obtained by summing up the two available traits
  values = sapply(multiPop@pops, function(pop){
    mean(pop@miscPop$pheno[,1:2])
  })
  expect_identical(
    selectPop(multiPop, nPop = 1, simParam = SP, trait = 1:2,
              use = function(pop, trait = trait) mean(pop@miscPop$pheno[,trait]))@pops[[1]]@id,
    multiPop@pops[[which.max(values)]]@id
  )
  expect_identical(
    selectPop(multiPop, nPop = 1, selectTop = F, simParam = SP, trait = 1:2,
              use = function(pop, trait = trait) mean(pop@miscPop$pheno[,trait]))@pops[[1]]@id,
    multiPop@pops[[which.min(values)]]@id
  )
  
  # Use a custom trait and use functions
  values = lapply(multiPop@pops, function(pop, trait = 1:2){
    c(pop@miscPop$pheno[,trait])
  })
  values = do.call('rbind', values)
  response = apply(values, 1, sd)
  
  expect_identical(
    selectPop(multiPop, nPop = 1, simParam = SP, selectTop = F,
              use = function(pop, trait = 1:2) c(pop@miscPop$pheno[,trait]),
              trait = function(x) apply(x, 1, sd))@pops[[1]]@id,
    multiPop@pops[[which.min(response)]]@id
  )
  expect_identical(
    selectPop(pop = multiPop, nPop = 1, simParam = SP, selectTop = T,
              use = function(pop, trait = 1:2) c(pop@miscPop$pheno[,trait]),
              trait = function(x) apply(x, 1, sd))@pops[[1]]@id,
    multiPop@pops[[which.max(response)]]@id
  )
  
  # MultiPop with 1 nested object
  
  multiPop = newMultiPop(pop[1:20], pop[21:40],
                         newMultiPop(pop[41:60],
                                     newMultiPop(pop[61:80], pop[81:100])))
  
  # setPhenoPop for all traits
  multiPop = setPhenoPop(multiPop, simParam = SP)
  
  # Selection can only be performed when all populations are Pop-class
  expect_error(
    selectPop(multiPop, nPop = 1, simParam = SP),
    paste(
      "This level contains", sum(sapply(multiPop@pops, isMultiPop)), 
      " MultiPop-class objects.",
      "\nSelection can only be performed when all populations at this level are",
      "Pop-class objects.\nYou may want to increase the value of 'level'"
    ), fixed = TRUE)
  expect_error(
    selectPop(multiPop, nPop = 1, level = 2, simParam = SP),
    paste(
      "This level contains", sum(sapply(multiPop@pops, isMultiPop)), 
      " MultiPop-class objects.",
      "\nSelection can only be performed when all populations at this level are",
      "Pop-class objects.\nYou may want to increase the value of 'level'"
    ), fixed = TRUE)
  
  # MultiPop with >1 nested object
  
  multiPop = newMultiPop(pop[1:20],
                         newMultiPop(pop[21:30],
                                     newMultiPop(pop[31:35], pop[36:40])),
                         newMultiPop(pop[41:60],
                                     newMultiPop(pop[61:80], pop[81:100])))
  
  multiPop = setPhenoPop(multiPop, simParam = SP)
  
  # Selection can only be performed when all populations are Pop-class
  expect_error(
    selectPop(multiPop, nPop = 1, simParam = SP),
    paste(
      "This level contains", sum(sapply(multiPop@pops, isMultiPop)), 
      " MultiPop-class objects.",
      "\nSelection can only be performed when all populations at this level are",
      "Pop-class objects.\nYou may want to increase the value of 'level'"
    ), fixed = TRUE)
  expect_error(
    selectPop(multiPop, nPop = 1, level = 2, simParam = SP),
    paste(
      "This level contains", sum(sapply(multiPop@pops[[3]]@pops, isMultiPop)), 
      " MultiPop-class objects.",
      "\nSelection can only be performed when all populations at this level are",
      "Pop-class objects.\nYou may want to increase the value of 'level'"
    ), fixed = TRUE)
})
