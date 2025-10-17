context("phenotypes")

test_that("asCategorical_converts_correctly",{
  cont = matrix(data = 0, nrow = 7, ncol = 3)
  cont[, 1] = c(-3, -2, -1, 0, 1, 2, 3)
  cont[, 2] = c(-3, -2, -1, 0, 1, 2, 3)
  cont[, 3] = c(-3, -2, -1, 0, 1, 2, 3)

  expect_equal(asCategorical(x = cont[, 1]),
               matrix(c(1, 1, 1, 2, 2, 2, 2)))
  expect_equal(asCategorical(x = cont[, 1], threshold = c(-1, 0, 1)),
               matrix(c(NA, NA, 1, 2, 2, NA, NA)))
  expect_equal(asCategorical(x = cont[, 1], threshold = c(-Inf, -1, 0, 1, Inf)),
               matrix(c(1, 1, 2, 3, 4, 4, 4)))

  expect_warning(asCategorical(x = cont[, 1], p = 0.5))
  expect_equal(suppressWarnings(asCategorical(x = cont[, 1], p = 0.5)),
               asCategorical(x = cont[, 1], p = c(0.5, 0.5)))

  trtMean = apply(X = cont, MARGIN = 2, FUN = mean)
  trtVar = apply(X = cont, MARGIN = 2, FUN = var)
  expect_equal(asCategorical(x = cont[, 1], p = c(0.5, 0.5), var = trtVar[1]),
               matrix(c(1, 1, 1, 2, 2, 2, 2)))
  expect_equal(asCategorical(x = cont[, 1], p = c(2/7, 1/7, 1/7, 3/7), var = trtVar[1]),
               matrix(c(1, 1, 2, 3, 4, 4, 4)))

  expect_error(asCategorical(x = cont))
  cont2 = asCategorical(x = cont,
                        threshold = list(NULL,
                                         c(-Inf, 0, Inf),
                                         NULL))
  cont2Exp = cont
  cont2Exp[, 2] = c(1, 1, 1, 2, 2, 2, 2)
  expect_equal(cont2, cont2Exp)

  expect_error(asCategorical(x = cont, p = c(0.5, 0.5)))
  pList = list(NULL, c(0.5, 0.5), NULL)
  expect_error(asCategorical(x = cont, p = pList))
  expect_error(asCategorical(x = cont, p = pList, mean = trtMean))
  cont2 = asCategorical(x = cont, p = pList, mean = trtMean, var=trtVar)
  cont2Exp = cont
  cont2Exp[, 2] = c(1, 1, 1, 2, 2, 2, 2)
  expect_equal(cont2, cont2Exp)
})

test_that("setPhenoPop", {
  getPhenoPop <- function(x) x@miscPop$pheno
  
  # Create founder haplotypes
  founderPop <- quickHaplo(nInd = 4, nChr = 1, segSites = 10)
  
  # Set simulation parameters and a single additive trait
  SP <- SimParam$new(founderPop)
  SP$addTraitA(10, mean = c(0, 0), var = c(1, 1))
  
  # Create a population and a multi-population object
  pop <- newPop(founderPop, simParam = SP)
  multiPop <- newMultiPop(pop[1:2], pop[3:4])
  
  # 1) Missing individual phenotypes should error when force = FALSE
  expect_error(
    setPhenoPop(multiPop, force = FALSE, simParam = SP),
    "The phenotypic matrix is empty. Use force=TRUE to create it."
  )
  
  # Create individual phenotypes (error variance = 1)
  multiPop <- setPheno(multiPop, varE = c(1, 1), simParam = SP)
  
  # After newMultiPop the miscPop slot should initially be empty
  expect_equal(length(multiPop@pops[[1]]@miscPop), 0)
  
  # 2) Set only trait 1 -> trait 2 columns should remain NA in miscPop
  multiPop <- setPhenoPop(multiPop, force = FALSE, traits = 1, simParam = SP)
  # Use onlyPhenoPop to get a clean matrix: one row per population, columns = traits
  phenoMat1 <- lapply(multiPop@pops, function(x) x@miscPop$pheno)
  phenoMat1 <- do.call('rbind', phenoMat1)
  expect_equal(dim(phenoMat1), c(length(multiPop@pops), SP$nTraits))
  # trait 2 should still be NA for all populations
  expect_true(all(is.na(phenoMat1[, 2])))
  
  # 3) Now set trait 2 using a different summary function (colSums),
  #    ensure trait 1 values are preserved
  multiPop <- setPhenoPop(multiPop, FUN = colSums, force = FALSE, traits = 2, simParam = SP)
  phenoMat2 <- lapply(multiPop@pops, function(x) x@miscPop$pheno)
  phenoMat2 <- do.call('rbind', phenoMat2)
  # trait 1 unchanged (still equal to previous phenoMat1 trait1)
  expect_equal(phenoMat2[, 1], phenoMat1[, 1])
  # trait 2 is now not all NA
  expect_false(all(is.na(phenoMat2[, 2])))
  
  # 4) Compare results when individual phenotypes are (re)generated via force = TRUE
  phenoMat_no_force <- setPhenoPop(multiPop, onlyPhenoPop = TRUE, simParam = SP)
  phenoMat_force <- setPhenoPop(multiPop, force = TRUE, varE = c(1, 1), onlyPhenoPop = TRUE, simParam = SP)
  expect_equal(dim(phenoMat_force), c(length(multiPop@pops), SP$nTraits))
  # In general the forced generation produces different values (random)
  expect_false(isTRUE(all.equal(phenoMat_no_force, phenoMat_force)))
  
  # 5) Custom function that returns medians; verify equality with direct medians
  med_fun <- function(x) apply(x, 2, median)
  phenoMed <- setPhenoPop(multiPop, onlyPhenoPop = TRUE, FUN = med_fun, simParam = SP)
  expect_equal(nrow(phenoMed), length(multiPop@pops))
  expect_equal(ncol(phenoMed), SP$nTraits)
  # compare for each population and trait
  for (i in seq_along(multiPop@pops)) {
    for (t in seq_len(SP$nTraits)) {
      expect_equal(unname(phenoMed[i,t]), median(multiPop@pops[[i]]@pheno[,t]))
    }
  }
  
  # 6) FUN returning a plain numeric vector (length == nTraits) should work
  vec_fun <- function(ph) { as.numeric(colMeans(ph)) }  # returns numeric vector
  phenoVec <- setPhenoPop(multiPop, onlyPhenoPop = TRUE, FUN = vec_fun, simParam = SP)
  expect_equal(dim(phenoVec), c(length(multiPop@pops), SP$nTraits))
  
  # 7) Calling with onlyPhenoPop on single Pop returns 1 x nTraits
  single_pheno <- setPhenoPop(multiPop@pops[[1]], onlyPhenoPop = TRUE, FUN = colMeans, simParam = SP)
  expect_equal(dim(single_pheno), c(1, SP$nTraits))
  
  # 8) Error handling: bad FUN output (wrong length) should raise an informative error
  bad_fun <- function(x) rep(1, SP$nTraits + 1)
  expect_error(setPhenoPop(multiPop, FUN = bad_fun, simParam = SP),
               "FUN returned an object of unexpected dimensions")
  
  # 9) Invalid input type should error
  expect_error(setPhenoPop(123, simParam = SP), 
               "pop must be an object of Pop, HybridPop or MultiPop class")
})
