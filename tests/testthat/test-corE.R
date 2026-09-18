context("corE")

# setPheno's corE argument builds a correlated error structure across traits.
# Nothing in the suite passed corE, so neither the correlation itself nor the
# checks guarding it were ever run.
#
# Correlation is a property of a sample, so most of these use a population
# large enough that the sample correlation lands close to what was asked for.
# The tolerances are loose on purpose: a wrong sign or a transposed matrix
# fails them, sampling noise does not.

corPop = function(nInd=2000, nChr=2, segSites=20, nTraits=2, corA=NULL,
                  seed=9951){
  set.seed(seed)
  founderPop = quickHaplo(nInd=nInd, nChr=nChr, segSites=segSites)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  if(is.null(corA)){
    corA = diag(nTraits)
  }
  SP$addTraitA(nQtlPerChr=10, mean=rep(0, nTraits), var=rep(1, nTraits),
               corA=corA)
  SP$setVarE(h2=rep(0.5, nTraits))
  pop = newPop(founderPop, simParam=SP)
  return(list(pop=pop, SP=SP))
}

# The correlation of the residuals, which is what corE controls
residCor = function(pop){
  resid = pheno(pop) - gv(pop)
  return(cor(resid[,1], resid[,2]))
}

test_that("corE puts the requested correlation into the residuals", {
  d = corPop()

  for(target in c(-0.6, 0, 0.6)){
    corE = matrix(c(1, target, target, 1), nrow=2)
    phenotyped = setPheno(d$pop, varE=c(1,1), corE=corE, simParam=d$SP)
    expect_equal(residCor(phenotyped), target, tolerance=0.08,
                 info=paste("target", target))
  }
})

test_that("corE leaves the error variances alone", {
  # Sampling based, and covered in outline by the block above
  skip_on_cran()
  d = corPop()
  corE = matrix(c(1, 0.5, 0.5, 1), nrow=2)
  phenotyped = setPheno(d$pop, varE=c(4, 9), corE=corE, simParam=d$SP)

  resid = pheno(phenotyped) - gv(phenotyped)
  expect_equal(var(resid[,1]), 4, tolerance=0.15)
  expect_equal(var(resid[,2]), 9, tolerance=0.15)

  # The correlation is unchanged by the variances being unequal
  expect_equal(cor(resid[,1], resid[,2]), 0.5, tolerance=0.08)
})

test_that("corE works from a variance covariance matrix as well as a vector", {
  # Sampling based
  skip_on_cran()
  d = corPop()
  corE = matrix(c(1, 0.4, 0.4, 1), nrow=2)

  # varE given as a matrix has its diagonal taken before corE is applied
  fromMatrix = setPheno(d$pop, varE=diag(c(4, 9)), corE=corE, simParam=d$SP)
  resid = pheno(fromMatrix) - gv(fromMatrix)

  expect_equal(var(resid[,1]), 4, tolerance=0.15)
  expect_equal(var(resid[,2]), 9, tolerance=0.15)
  expect_equal(cor(resid[,1], resid[,2]), 0.4, tolerance=0.08)
})

test_that("correlated errors do not disturb the genetic values", {
  # Sampling based
  skip_on_cran()
  d = corPop()
  corE = matrix(c(1, 0.7, 0.7, 1), nrow=2)

  before = gv(d$pop)
  phenotyped = setPheno(d$pop, varE=c(1,1), corE=corE, simParam=d$SP)
  expect_equal(gv(phenotyped), before)

  # And the phenotypic mean is still the genetic mean
  expect_equal(colMeans(pheno(phenotyped)), colMeans(before), tolerance=0.1)
})

test_that("corE is checked for shape and symmetry", {
  d = corPop(nInd=50)

  # Not symmetric
  expect_error(setPheno(d$pop, varE=c(1,1),
                        corE=matrix(c(1, 0.5, 0.2, 1), nrow=2),
                        simParam=d$SP),
               "symmetric")

  # Wrong size for the number of traits
  expect_error(setPheno(d$pop, varE=c(1,1), corE=diag(3), simParam=d$SP),
               "corE")
})

test_that("varE is checked against the number of traits", {
  d = corPop(nInd=50)

  expect_error(setPheno(d$pop, varE=c(1,1,1), simParam=d$SP), "varE")
  expect_error(setPheno(d$pop, varE=matrix(1, nrow=3, ncol=3),
                        simParam=d$SP), "varE")
  expect_error(setPheno(d$pop, varE=matrix(c(1,0.5,0.2,1), nrow=2),
                        simParam=d$SP), "symmetric")
})

test_that("a genetic correlation and an error correlation are separate", {
  # Sampling based, and builds a second population of its own
  skip_on_cran()
  # Genetically correlated traits with uncorrelated errors, and the reverse.
  # Each should show up in its own place and not the other.
  corA = matrix(c(1, 0.8, 0.8, 1), nrow=2)
  d = corPop(corA=corA, seed=9952)

  genCor = cor(gv(d$pop)[,1], gv(d$pop)[,2])
  expect_gt(genCor, 0.5)

  # Uncorrelated errors leave the residual correlation near zero
  indep = setPheno(d$pop, varE=c(1,1), corE=diag(2), simParam=d$SP)
  expect_equal(residCor(indep), 0, tolerance=0.08)

  # Negatively correlated errors show up in the residuals without changing
  # the genetic correlation
  opposed = setPheno(d$pop, varE=c(1,1),
                     corE=matrix(c(1,-0.7,-0.7,1), nrow=2), simParam=d$SP)
  expect_equal(residCor(opposed), -0.7, tolerance=0.08)
  expect_equal(cor(gv(opposed)[,1], gv(opposed)[,2]), genCor)
})

test_that("multi-trait setPheno handles reps and h2 together", {
  # Sampling based, with four phenotyping passes over the population
  skip_on_cran()
  d = corPop()

  # Repeated measurement divides the error variance by the number of reps
  once = setPheno(d$pop, varE=c(4,4), reps=1, simParam=d$SP)
  four = setPheno(d$pop, varE=c(4,4), reps=4, simParam=d$SP)

  varOnce = var(pheno(once)[,1] - gv(once)[,1])
  varFour = var(pheno(four)[,1] - gv(four)[,1])
  expect_equal(varOnce, 4, tolerance=0.15)
  expect_equal(varFour, 1, tolerance=0.15)

  # h2 given per trait sets the error variance to match
  byH2 = setPheno(d$pop, h2=c(0.25, 0.75), simParam=d$SP)
  resid = pheno(byH2) - gv(byH2)
  h2Obs = apply(gv(byH2), 2, var) /
    (apply(gv(byH2), 2, var) + apply(resid, 2, var))
  expect_equal(unname(h2Obs[1]), 0.25, tolerance=0.12)
  expect_equal(unname(h2Obs[2]), 0.75, tolerance=0.12)
})

test_that("onlyPheno returns the phenotypes instead of a population", {
  d = corPop(nInd=100)
  out = setPheno(d$pop, varE=c(1,1), corE=diag(2), onlyPheno=TRUE,
                 simParam=d$SP)

  expect_true(is.matrix(out))
  expect_false(isPop(out))
  expect_equal(nrow(out), nInd(d$pop))
  expect_equal(ncol(out), 2L)
  expect_true(all(is.finite(out)))
})
